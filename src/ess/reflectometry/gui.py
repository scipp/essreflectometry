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

    def run_table_from_metadata(self, df):
        raise NotImplementedError()

    def reduction_table_from_run_table(self, df):
        raise NotImplementedError()

    def reference_table_from_run_table(self, df):
        raise NotImplementedError()

    def display_results(self, results):
        raise NotImplementedError()

    def run_workflow(self):
        raise NotImplementedError()

    def log(self, message):
        with self.output:
            display(message)

    @property
    def path(self):
        p = self.proposal_number_box.value
        if p.isdigit():
            # Handling proposal numbers is not yet implemented
            raise NotImplementedError()
        return p

    def __init__(self):
        self.output = widgets.Output()
        self.log("init")

        self.runs_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            layout={"height": "300px"},
            column_visibility={"key": False},
        )
        self.reduction_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            layout={"height": "300px"},
            column_visibility={"key": False},
        )
        self.reference_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            layout={"height": "300px"},
            column_visibility={"key": False},
        )

        header1 = widgets.Label("Runs")
        header2 = widgets.Label("Reduction")
        header3 = widgets.Label("Reference runs")

        self.proposal_number_box = widgets.Text(
            value="",
            placeholder="Proposal number or file path",
            description="Proposal no.:",
            layout=widgets.Layout(description_width="auto"),
            disabled=False,
        )

        load_files_button = widgets.Button(description="Load runs")
        transfer_button = widgets.Button(description="Transfer")
        reduce_button = widgets.Button(description="Reduce")

        def load_files(b):
            self.log("load files")
            self.runs_table.data = self.run_table_from_metadata(
                self.load_runs_from_path(self.path)
            )

        def transfer_runs(b):
            self.log("transfer runs")
            self.reference_table.data = self.reference_table_from_run_table(
                self.runs_table.data
            )
            self.reduction_table.data = self.reduction_table_from_run_table(
                self.runs_table.data
            )

        def reduce_data(b):
            self.log("reduce data")
            self.display_results(self.run_workflow())

        load_files_button.on_click(load_files)
        transfer_button.on_click(transfer_runs)
        reduce_button.on_click(reduce_data)

        data_buttons = widgets.HBox([load_files_button, transfer_button, reduce_button])

        tab_data = widgets.VBox(
            [
                data_buttons,
                widgets.HBox(
                    [
                        widgets.VBox(
                            [
                                header1,
                                self.runs_table,
                            ],
                            layout={"width": "100%"},
                        ),
                        widgets.VBox(
                            [
                                header2,
                                self.reduction_table,
                            ],
                            layout={"width": "100%"},
                        ),
                    ]
                ),
            ],
        )
        tab_settings = widgets.VBox(
            [widgets.Label("This is the settings tab"), header3, self.reference_table],
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
                self.output,
            ]
        )

    def load_runs_from_path(self, path):
        self.log("load runs from path")
        metadata = [
            self.read_meta_data(fpath)
            for fpath in glob.glob(os.path.join(path, '*.hdf'))
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

    def run_table_from_metadata(self, df):
        return df.sort_values(by='Run')

    def reduction_table_from_run_table(self, df):
        return (
            df[df["Sample"] != "sm5"]
            .groupby(["Sample", "Angle"], as_index=False)
            .agg(Runs=("Run", list))
            .sort_values(by=["Sample", "Angle"])
        )

    def reference_table_from_run_table(self, df):
        return (
            df[df["Sample"] == "sm5"]
            .groupby(["Sample"], as_index=False)
            .agg(Runs=("Run", list))
            .sort_values(by="Sample")
        )

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
        with self.output:
            display(tiled)

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

        workflow[QBins] = sc.geomspace(
            dim='Q', start=0.005, stop=0.4, num=391, unit='1/angstrom'
        )
        workflow[WavelengthBins] = sc.geomspace(
            'wavelength', 3, 12.5, 2001, unit='angstrom'
        )

        workflow[YIndexLimits] = sc.scalar(17), sc.scalar(47)
        workflow[ZIndexLimits] = sc.scalar(60), sc.scalar(380)

        progress = widgets.IntProgress(min=0, max=len(sample_df))

        with self.output:
            display(progress)

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
            reflectivity_curves.append(wf.compute(ReflectivityOverQ).hist())
            progress.value += 1

        return reflectivity_curves
