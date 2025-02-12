import glob
import os

import h5py
import ipywidgets as widgets
import pandas as pd
import plopp as pp
import scipp as sc
from ipydatagrid import DataGrid
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

    def display_results(self, results):
        raise NotImplementedError()

    def run_workflow(self):
        raise NotImplementedError()

    def log(self, message):
        out = widgets.Output()
        with out:
            display(message)
        self.logbox.children = (out, *self.logbox.children)

    def sync(self, *_):
        db = {}
        # db["settings"] = self.load_settings()
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

    @property
    def path(self):
        if self._path is None:
            raise ValueError("Path is not set")
        return self._path

    def __init__(self):
        self.logbox = widgets.VBox([])
        self._path = None
        self.log("init")

        self.runs_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
        )
        self.reduction_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
        )
        self.reference_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
        )
        self.runs_table.on_cell_change(self.sync)
        self.reduction_table.on_cell_change(self.sync)
        self.reference_table.on_cell_change(self.sync)

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

        def reduce_data(_):
            self.log("reduce data")
            self.display_results(self.run_workflow())

        reduce_button.on_click(reduce_data)

        add_row_button = widgets.Button(description="Add row")
        delete_row_button = widgets.Button(description="Remove row")

        def add_row(_):
            self.log("add row")
            row = self.runs_table.data.iloc[-1:].copy()
            row[row.columns[1]] = row[row.columns[1]] + '_copy'
            self.runs_table.data = pd.concat([self.runs_table.data, row])
            self.sync()

        def delete_row(_):
            self.log("delete row")
            self.runs_table.data = self.runs_table.data.iloc[:-1]
            self.sync()

        add_row_button.on_click(add_row)
        delete_row_button.on_click(delete_row)

        data_buttons = widgets.HBox([add_row_button, delete_row_button, reduce_button])

        tab_data = widgets.VBox(
            [
                data_buttons,
                widgets.HBox(
                    [
                        widgets.VBox(
                            [
                                widgets.Label("Runs"),
                                self.runs_table,
                            ],
                            layout={"width": "100%"},
                        ),
                        widgets.VBox(
                            [
                                widgets.Label("Reduction"),
                                self.reduction_table,
                            ],
                            layout={"width": "100%"},
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
    def _merge_old_and_new_state(new, old, on):
        old = old if on in old else old.assign(**{on: None})
        new = new if on in new else new.assign(**{on: None})
        df = new.merge(old, how='left', on=on)
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
        user_reduction = db['user_reduction'].drop(
            columns=["Sample", "Angle"], errors='ignore'
        )
        df = self._merge_old_and_new_state(df, user_reduction, on='Runs')
        df = df.drop_duplicates(("Sample", "Angle", "Runs"))
        self._setdefault(df, "QBins", 391)
        self._setdefault(df, "QStart", 0.01)
        self._setdefault(df, "QStop", 0.3)
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
        user_reference = db['user_reference'].drop(columns=["Sample"], errors='ignore')
        df = self._merge_old_and_new_state(df, user_reference, on='Runs')
        df = self._ordercolumns(df, 'Sample')
        return df.sort_values(by="Sample")

    def display_results(self, results):
        df = self.reduction_table.data.copy()
        df["rownum"] = range(len(df))
        to_combine = df.groupby("Sample", as_index=False).agg({"rownum": list})
        tiled = pp.tiled(1, 2)
        tiled[0, 0] = pp.plot(
            {
                ','.join(params["Runs"]): res
                for (_, params), res in zip(df.iterrows(), results, strict=True)
            },
            norm='log',
            figsize=(10, 7),
        )
        tiled[0, 1] = pp.plot(
            {
                params["Sample"]: combine_curves(
                    [results[i] for i in params['rownum']],
                    q_bin_edges=results[0].coords['Q'],
                )
                for _, params in to_combine.iterrows()
            },
            norm='log',
            figsize=(10, 7),
        )
        self.log(tiled)

    def get_filepath_from_run(self, run):
        return os.path.join(self.path, f'amor2024n{run:0>6}.hdf')

    def run_workflow(self):
        sample_df = self.reduction_table.data
        reference_df = self.reference_table.data

        workflow = amor.AmorWorkflow()
        workflow[SampleSize[SampleRun]] = sc.scalar(10, unit='mm')
        workflow[SampleSize[ReferenceRun]] = sc.scalar(10, unit='mm')

        workflow[ChopperPhase[ReferenceRun]] = sc.scalar(7.5, unit='deg')
        workflow[ChopperPhase[SampleRun]] = sc.scalar(7.5, unit='deg')

        workflow[WavelengthBins] = sc.geomspace(
            'wavelength', 3, 12.5, 2001, unit='angstrom'
        )

        workflow[YIndexLimits] = sc.scalar(17), sc.scalar(47)
        workflow[ZIndexLimits] = sc.scalar(60), sc.scalar(380)

        progress = widgets.IntProgress(min=0, max=len(sample_df))
        self.log(progress)

        runs = (
            reference_df.iloc[0]["Runs"]
            if not isinstance(reference_df.iloc[0]["Runs"], str)
            else reference_df.iloc[0]["Runs"].split(',')
        )
        workflow[ReducedReference] = with_filenames(
            workflow, ReferenceRun, list(map(self.get_filepath_from_run, runs))
        ).compute(ReducedReference)
        progress.value += 1

        reflectivity_curves = []
        for _, params in sample_df.iterrows():
            runs = (
                params["Runs"]
                if not isinstance(params['Runs'], str)
                else params['Runs'].split(',')
            )
            wf = with_filenames(
                workflow, SampleRun, list(map(self.get_filepath_from_run, runs))
            )
            wf[QBins] = sc.geomspace(
                dim='Q',
                start=params['QStart'],
                stop=params['QStop'],
                num=int(params['QBins']),
                unit='1/angstrom',
            )
            reflectivity_curves.append(wf.compute(ReflectivityOverQ).hist())
            progress.value += 1

        return reflectivity_curves
