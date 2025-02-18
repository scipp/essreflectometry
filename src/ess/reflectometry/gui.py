import glob
import os

import h5py
import ipywidgets as widgets
import matplotlib
import pandas as pd
import plopp as pp
import scipp as sc
from ipydatagrid import DataGrid, VegaExpr
from IPython.display import display

from ess import amor
from ess.amor.types import ChopperPhase
from ess.reflectometry.tools import combine_curves
from ess.reflectometry.types import (
    QBins,
    ReducedReference,
    ReferenceRun,
    ReflectivityOverQ,
    SampleRun,
    SampleSize,
    WavelengthBins,
    YIndexLimits,
    ZIndexLimits,
)
from ess.reflectometry.workflow import with_filenames


class ReflectometryBatchReductionGUI:
    def read_meta_data(self, path):
        raise NotImplementedError()

    def sync_runs_table(self, db):
        raise NotImplementedError()

    def sync_reduction_table(self, db):
        raise NotImplementedError()

    def sync_reference_table(self, db):
        raise NotImplementedError()

    def sync_custom_reduction_table(self):
        raise NotImplementedError()

    def display_results(self):
        raise NotImplementedError()

    def run_workflow(self):
        raise NotImplementedError()

    def get_row_key(self, row):
        reference_metadata = (
            tuple(self.reference_table.data.iloc[0])
            if len(self.reference_table.data) > 0
            else (None,)
        )
        return (tuple(row), tuple(reference_metadata))

    def sync_table_colors(self, table):
        template = 'row == {i} ? {reduced_color} : '
        expr = ''
        for i, (_, row) in enumerate(table.data.iterrows()):
            for row_key in self.results.keys():
                if self.get_row_key(row) == row_key:
                    expr += template.format(i=i, reduced_color="'lightgreen'")
        expr += "default_value"
        table.default_renderer.background_color = VegaExpr(expr)

    def set_result(self, metadata, result):
        self.results[self.get_row_key(metadata)] = result
        self.sync_table_colors(self.reduction_table)
        self.sync_table_colors(self.custom_reduction_table)
        self.sync_table_colors(self.reference_table)

    def log(self, message):
        out = widgets.Output()
        with out:
            display(message)
        self.logbox.children = (out, *self.logbox.children)

    @staticmethod
    def set_table_height(table, extra=0):
        height = (len(table.data) + 1) * (table.base_row_size + 1) + 5 + extra
        table.layout.height = f'{height}px'

    def sync(self, *_):
        db = {}
        # db["settings"] = self.load_settings()
        db["run_number_min"] = int(self.run_number_min.value)
        db["run_number_max"] = int(self.run_number_max.value)
        db["meta"] = self.load_runs()
        db["user_runs"] = self.runs_table.data
        db["user_reduction"] = self.reduction_table.data
        db["user_reference"] = self.reference_table.data

        db["user_runs"] = self.sync_runs_table(db)
        db["user_reduction"] = self.sync_reduction_table(db)
        db["user_reference"] = self.sync_reference_table(db)

        self.runs_table.data = db["user_runs"]
        self.reduction_table.data = db["user_reduction"]
        self.reference_table.data = db["user_reference"]

        self.set_table_height(self.runs_table)
        self.set_table_height(self.reduction_table)
        self.set_table_height(self.custom_reduction_table)
        self.sync_table_colors(self.reduction_table)
        self.sync_table_colors(self.custom_reduction_table)
        self.sync_table_colors(self.reference_table)

    @property
    def path(self):
        if self._path is None:
            raise ValueError("Path is not set")
        return self._path

    def __init__(self):
        self.logbox = widgets.VBox([])
        self._path = None
        self.log("init")

        self.results = {}

        self.runs_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
            selection_mode="row",
        )
        self.reduction_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
            selection_mode="row",
        )
        self.reference_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
            selection_mode="row",
        )
        self.custom_reduction_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
            selection_mode="row",
        )

        self.runs_table.on_cell_change(self.sync)
        self.reduction_table.on_cell_change(self.sync)
        self.reference_table.on_cell_change(self.sync)

        self.custom_reduction_table.on_cell_change(
            lambda _: self.sync_custom_reduction_table()
        )

        self.proposal_number_box = widgets.Text(
            value="",
            placeholder="Proposal number or file path",
            description="Proposal no.:",
            layout=widgets.Layout(description_width="auto"),
            disabled=False,
        )

        def set_proposal_number_state(state):
            if state == "good":
                self.proposal_number_box.layout.border = '2px solid green'
            if state == "bad":
                self.proposal_number_box.layout.border = '2px solid red'

        def on_proposal_number_change(_):
            p = self.proposal_number_box.value
            if p.isdigit():
                # Handling proposal numbers is not yet implemented
                self._path = None
                set_proposal_number_state("bad")
            elif not os.path.isdir(p):
                self._path = None
                set_proposal_number_state("bad")
            else:
                self._path = p
                set_proposal_number_state("good")
                self.sync()

        set_proposal_number_state("bad")
        self.proposal_number_box.observe(on_proposal_number_change, names='value')

        reduce_button = widgets.Button(description="Reduce")
        plot_button = widgets.Button(description="Plot")

        def reduce_data(_):
            self.log("reduce data")
            self.run_workflow()

        def plot_results(_):
            self.log("plot results")
            self.display_results()

        reduce_button.on_click(reduce_data)
        plot_button.on_click(plot_results)

        add_row_button = widgets.Button(description="Add row")
        delete_row_button = widgets.Button(description="Remove row")

        def add_row(_):
            self.log("add row")
            if len(self.custom_reduction_table.data) > 0:
                row = self.custom_reduction_table.data.iloc[-1:]
            elif len(self.reduction_table.data) > 0:
                row = self.reduction_table.data.iloc[-1:]
            else:
                return
            # To avoid a flickering scrollbar
            self.set_table_height(self.custom_reduction_table, extra=25)
            self.custom_reduction_table.data = pd.concat(
                [self.custom_reduction_table.data, row]
            )
            self.set_table_height(self.custom_reduction_table)
            self.sync_table_colors(self.custom_reduction_table)

        def delete_row(_):
            self.log("delete row")
            self.custom_reduction_table.data = self.custom_reduction_table.data.iloc[
                :-1
            ]
            self.set_table_height(self.custom_reduction_table)

        add_row_button.on_click(add_row)
        delete_row_button.on_click(delete_row)
        data_buttons = widgets.HBox(
            [add_row_button, delete_row_button, reduce_button, plot_button]
        )

        self.run_number_min = widgets.IntText(
            value=0, description='', layout=widgets.Layout(width='5em')
        )
        self.run_number_max = widgets.IntText(
            value=9999, description='', layout=widgets.Layout(width='5em')
        )
        self.run_number_min.observe(self.sync, names='value')
        self.run_number_max.observe(self.sync, names='value')

        run_number_filter = widgets.HBox(
            [self.run_number_min, widgets.Label("<=Run<="), self.run_number_max]
        )

        tab_data = widgets.VBox(
            [
                widgets.HBox(
                    [
                        widgets.VBox(
                            [
                                widgets.Label("Runs"),
                                run_number_filter,
                                self.runs_table,
                            ],
                            layout={"width": "35%"},
                        ),
                        widgets.VBox(
                            [
                                widgets.Label("Reduction"),
                                data_buttons,
                                self.reduction_table,
                                self.custom_reduction_table,
                            ],
                            layout={"width": "60%"},
                        ),
                    ]
                ),
            ],
        )
        tab_settings = widgets.VBox(
            [
                widgets.Label("This is the settings tab"),
                widgets.Label("Reference runs"),
                self.reference_table,
            ],
            layout={"width": "100%"},
        )

        self.tabs = widgets.Tab()
        self.tabs.children = [tab_data, tab_settings]
        self.tabs.set_title(0, "Reduce")
        self.tabs.set_title(1, "Settings")

        self.main = widgets.VBox(
            [
                self.proposal_number_box,
                self.tabs,
                self.logbox,
            ]
        )

    def load_runs(self):
        self.log("load runs from path")
        metadata = [
            self.read_meta_data(fpath)
            for fpath in glob.glob(os.path.join(self.path, '*.hdf'))
        ]
        return pd.DataFrame(metadata)

    @property
    def widget(self):
        return self.main


class AmorBatchReductionGUI(ReflectometryBatchReductionGUI):
    def read_meta_data(self, path):
        with h5py.File(path) as f:
            return {
                "Sample": f['entry1']['sample']['name'][()][0].decode('utf8'),
                "Run": path[-8:-4],
                "Angle": f['entry1']['Amor']['master_parameters']['mu']['value'][0, 0],
            }

    @staticmethod
    def _merge_old_and_new_state(new, old, on, how='left'):
        old = old if on in old else old.assign(**{on: None})
        new = new if on in new else new.assign(**{on: None})
        df = new.merge(old, how=how, on=on)
        for right in df.columns:
            if right.endswith("_y"):
                new = right.removesuffix("_y")
                left = new + "_x"
                df[new] = df[right].combine_first(df[left])
                df = df.drop(columns=[left, right])
        return df

    @staticmethod
    def _setdefault(df, col, value):
        df[col] = value if col not in df.columns else df[col].fillna(value)

    @staticmethod
    def _ordercolumns(df, *cols):
        columns = [*cols, *sorted(set(df.columns) - {*cols})]
        return df[columns]

    def sync_runs_table(self, db):
        df = self._merge_old_and_new_state(db["meta"], db["user_runs"], on='Run')
        df = df[db['run_number_min'] <= df['Run'].astype(int)][
            db['run_number_max'] >= df['Run'].astype(int)
        ]
        self._setdefault(df, "Exclude", False)
        df = self._ordercolumns(df, 'Run', 'Sample')
        return df.sort_values(by='Run')

    def sync_reduction_table(self, db):
        df = db["user_runs"]
        df = (
            df[df["Sample"] != "sm5"][~df["Exclude"]]
            .groupby(["Sample", "Angle"], as_index=False)
            .agg(Runs=("Run", tuple))
            .sort_values(["Sample", "Angle"])
        )
        # We don't want changes to Sample or Angle made
        # in the user_reduction table to persist
        user_reduction = db['user_reduction'].drop(
            columns=["Sample", "Angle"], errors='ignore'
        )
        df = self._merge_old_and_new_state(df, user_reduction, on='Runs')
        self._setdefault(df, "QBins", 391)
        self._setdefault(df, "QStart", 0.01)
        self._setdefault(df, "QStop", 0.3)
        self._setdefault(df, "Scale", 1.0)
        df = self._ordercolumns(df, 'Sample', 'Angle', 'Runs')
        return df.sort_values(["Sample", "Angle"])

    def sync_reference_table(self, db):
        df = db["user_runs"]
        df = (
            df[df["Sample"] == "sm5"][~df["Exclude"]]
            .groupby(["Sample"], as_index=False)
            .agg(Runs=("Run", tuple))
            .sort_values(by="Sample")
        )
        # We don't want changes to Sample
        # in the user_reference table to persist
        user_reference = db['user_reference'].drop(columns=["Sample"], errors='ignore')
        df = self._merge_old_and_new_state(df, user_reference, on='Runs')
        self._setdefault(df, "Ymin", 17)
        self._setdefault(df, "Ymax", 47)
        self._setdefault(df, "Zmin", 60)
        self._setdefault(df, "Zmax", 380)
        self._setdefault(df, "Lmin", 3.0)
        self._setdefault(df, "Lmax", 12.5)
        df = self._ordercolumns(df, 'Sample', 'Runs')
        return df.sort_values(by="Sample")

    def sync_custom_reduction_table(self):
        df = self.custom_reduction_table.data.copy()
        if 'Runs' in df.columns:
            df['Runs'] = df['Runs'].map(
                lambda x: tuple(x)
                if isinstance(x, tuple | list)
                else tuple(x.split(','))
            )
        self.custom_reduction_table.data = df
        self.sync_table_colors(self.custom_reduction_table)

    def display_results(self):
        df = self.get_selected_rows()
        try:
            results = [
                next(v for (m, _), v in self.results.items() if m == key)
                for key in (tuple(row) for _, row in df.iterrows())
            ]
        except StopIteration:
            # No results were found for the selected row
            # It hasn't been computed yet, so compute it and try again.
            self.run_workflow()
            self.display_results()
            return

        df["rownum"] = range(len(df))
        to_combine = df.groupby("Sample", as_index=False).agg({"rownum": list})

        def get_unique_names(df):
            names = [','.join(params["Runs"]) for (_, params) in df.iterrows()]
            duplicated_name_counter = {}
            unique = []
            for i, name in enumerate(names):
                if name not in names[:i]:
                    unique.append(name)
                else:
                    duplicated_name_counter.setdefault(name, 0)
                    duplicated_name_counter[name] += 1
                    unique.append(f'{name}_{duplicated_name_counter[name]}')
            return unique

        all_runs_plot = pp.plot(
            dict(zip(get_unique_names(df), results, strict=True)),
            norm='log',
            vmin=max(1e-6, min(result.min().value for result in results)),
        )

        def get_q_bin_edges(rows):
            qmin = min(sc.min(results[i].coords['Q']) for i in rows)
            qmax = max(sc.max(results[i].coords['Q']) for i in rows)
            qnum = 2 * max(results[i].coords['Q'].size for i in rows)
            return sc.linspace('Q', qmin, qmax, qnum)

        stitched_plot = pp.plot(
            {
                params["Sample"]: combine_curves(
                    [results[i] for i in params['rownum']],
                    q_bin_edges=get_q_bin_edges(params['rownum']),
                )
                for _, params in to_combine.iterrows()
            },
            norm='log',
            vmin=max(1e-6, min(result.min().value for result in results)),
        )
        if matplotlib.get_backend().lower().startswith('qt'):
            all_runs_plot.show()
            stitched_plot.show()
        else:
            tiled = pp.tiled(1, 2)
            tiled[0, 0] = all_runs_plot
            tiled[0, 1] = stitched_plot
            self.log(tiled)

    def get_filepath_from_run(self, run):
        return os.path.join(self.path, f'amor2024n{run:0>6}.hdf')

    def get_selected_rows(self):
        chunks = [
            table.data.iloc[s['r1'] : s['r2'] + 1]
            for table in (self.reduction_table, self.custom_reduction_table)
            for s in table.selections
        ]
        if len(chunks) == 0:
            chunks = [self.reduction_table.data, self.custom_reduction_table.data]
        return pd.concat(chunks)

    def run_workflow(self):
        sample_df = self.get_selected_rows()
        reference_df = self.reference_table.data.iloc[0]

        workflow = amor.AmorWorkflow()
        workflow[SampleSize[SampleRun]] = sc.scalar(10, unit='mm')
        workflow[SampleSize[ReferenceRun]] = sc.scalar(10, unit='mm')

        workflow[ChopperPhase[ReferenceRun]] = sc.scalar(7.5, unit='deg')
        workflow[ChopperPhase[SampleRun]] = sc.scalar(7.5, unit='deg')

        workflow[WavelengthBins] = sc.geomspace(
            'wavelength',
            reference_df['Lmin'],
            reference_df['Lmax'],
            2001,
            unit='angstrom',
        )

        workflow[YIndexLimits] = (
            sc.scalar(reference_df['Ymin']),
            sc.scalar(reference_df['Ymax']),
        )
        workflow[ZIndexLimits] = (
            sc.scalar(reference_df['Zmin']),
            sc.scalar(reference_df['Zmax']),
        )

        progress = widgets.IntProgress(min=0, max=len(sample_df))
        self.log(progress)

        if (key := self.get_row_key(reference_df)) in self.results:
            reference_result = self.results[key]
        else:
            reference_result = with_filenames(
                workflow,
                ReferenceRun,
                list(map(self.get_filepath_from_run, reference_df["Runs"])),
            ).compute(ReducedReference)
            self.set_result(reference_df, reference_result)

        workflow[ReducedReference] = reference_result
        progress.value += 1

        for _, params in sample_df.iterrows():
            if (key := self.get_row_key(params)) in self.results:
                progress.value += 1
                continue

            wf = with_filenames(
                workflow,
                SampleRun,
                list(map(self.get_filepath_from_run, params['Runs'])),
            )
            wf[QBins] = sc.geomspace(
                dim='Q',
                start=params['QStart'],
                stop=params['QStop'],
                num=int(params['QBins']),
                unit='1/angstrom',
            )
            self.set_result(
                params, params["Scale"] * wf.compute(ReflectivityOverQ).hist()
            )
            progress.value += 1
