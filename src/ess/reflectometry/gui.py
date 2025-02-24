import glob
import os

import h5py
import ipywidgets as widgets
import matplotlib
import numpy as np
import pandas as pd
import plopp as pp
import scipp as sc
from ipydatagrid import DataGrid, VegaExpr
from IPython.display import display
from ipytree import Node, Tree

from ess import amor
from ess.amor.types import ChopperPhase
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
    """GUI for batch reduction of reflectometry data.

    Known limitations:
    1. Remove plot button behavior is inconsistent:
       - Removes the target plot with its controls
       - Previous plots disappear but their control buttons remain
       - Previous plots maintain interactivity despite attempted conversion to static
    2. Dataset toggle does not affect error bars as they are separate matplotlib artists
    3. Remove row button removes last row instead of selected row
    4. LogY toggle doesn't work due to workarounds for plopp's axis behavior:
       - Plopp's autoscale was flipping the y-axis orientation
       - We override multiple plopp/matplotlib methods to maintain correct orientation
       - This prevents the LogY toggle from working as it would interfere with our fixes

    These limitations are documented with FIXME comments in the relevant code sections.
    """

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
        self.text_log.children = (out, *self.text_log.children)

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
        self.text_log = widgets.VBox([])
        self.progress_log = widgets.VBox([])
        self.plot_log = widgets.VBox([])
        self.plot_counter = 0  # Add a counter to create unique IDs for plots
        self._path = None
        self.log("init")

        self.results = {}

        self.runs_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
            selection_mode="cell",
        )
        self.reduction_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
            selection_mode="cell",
        )
        self.reference_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
            selection_mode="cell",
        )
        self.custom_reduction_table = DataGrid(
            pd.DataFrame([]),
            editable=True,
            auto_fit_columns=True,
            column_visibility={"key": False},
            selection_mode="cell",
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
            # Check if there are any selections in the reduction table
            if len(self.reduction_table.selections) > 0:
                # Get the first selected row
                selection = self.reduction_table.selections[0]
                row = self.reduction_table.data.iloc[
                    selection['r1'] : selection['r2'] + 1
                ]
            else:
                # Create an empty row with default values
                row = pd.DataFrame(
                    [
                        {
                            'Sample': '',
                            'Angle': 0.0,
                            'Runs': (),
                            'QBins': 391,
                            'QStart': 0.01,
                            'QStop': 0.3,
                            'Scale': 1.0,
                        }
                    ]
                )
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
        data_buttons = widgets.HBox([reduce_button, plot_button])

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
                                widgets.Label("Runs Table"),
                                run_number_filter,
                                self.runs_table,
                            ],
                            layout={"width": "35%"},
                        ),
                        widgets.VBox(
                            [
                                data_buttons,
                                widgets.VBox(
                                    [
                                        widgets.Label("Auto Reduction Table"),
                                        self.reduction_table,
                                    ],
                                    layout={'margin': '10px 0'},
                                ),
                                widgets.VBox(
                                    [
                                        widgets.Label("Manual Reduction Table"),
                                        widgets.HBox(
                                            [add_row_button, delete_row_button],
                                            layout={'margin': '5px 0'},
                                        ),
                                        self.custom_reduction_table,
                                    ],
                                    layout={'margin': '10px 0'},
                                ),
                            ],
                            layout={"width": "60%"},
                        ),
                    ]
                ),
                widgets.VBox(
                    [
                        widgets.VBox(
                            [widgets.Label("Progress"), self.progress_log],
                            layout={'width': '100%', 'margin': '10px 0'},
                        ),
                        widgets.VBox(
                            [widgets.Label("Plots"), self.plot_log],
                            layout={'width': '100%', 'margin': '10px 0'},
                        ),
                    ]
                ),
            ]
        )

        tab_settings = widgets.VBox(
            [
                widgets.Label("This is the settings tab"),
                widgets.Label("Reference runs"),
                self.reference_table,
            ],
            layout={"width": "100%"},
        )

        tab_log = widgets.VBox(
            [widgets.Label("Messages"), self.text_log],
            layout={"width": "100%"},
        )

        # Create tree widget for Nexus structure
        self.nexus_tree = Tree(
            layout=widgets.Layout(
                width='100%',
                height='100%',  # Fill the container height
            )
        )
        self.nexus_tree.nodes = [Node("Select a run to view its structure")]

        # Add selection handler to runs table
        self.runs_table.observe(self.update_nexus_view, names='selections')

        # Create content viewer widget
        self.nexus_content = widgets.Textarea(
            value='Select a node to view its content',
            layout=widgets.Layout(width='100%', height='600px'),
            disabled=True,  # Make it read-only
        )

        # Add selection handler to tree
        self.nexus_tree.observe(self.on_tree_select, names='selected_nodes')

        # Create the Nexus Explorer tab content
        tab_nexus = widgets.VBox(
            [
                widgets.Label("Nexus Explorer"),
                widgets.HBox(
                    [
                        widgets.VBox(
                            [
                                widgets.Label("Runs Table"),
                                self.runs_table,
                            ],
                            layout={"width": "30%"},
                        ),
                        widgets.VBox(
                            [
                                widgets.Label("File Structure"),
                                widgets.VBox(
                                    [self.nexus_tree],
                                    layout=widgets.Layout(
                                        width='100%',
                                        height='600px',
                                        min_height='100px',  # Min resize height
                                        max_height='1000px',  # Max resize height
                                        overflow_y='scroll',
                                        border='1px solid lightgray',
                                        resize='vertical',  # Add resize handle
                                    ),
                                ),
                            ],
                            layout={"width": "35%"},
                        ),
                        widgets.VBox(
                            [
                                widgets.Label("Content"),
                                self.nexus_content,
                            ],
                            layout={"width": "35%"},
                        ),
                    ]
                ),
            ],
            layout={"width": "100%"},
        )

        self.tabs = widgets.Tab()
        self.tabs.children = [
            tab_data,
            tab_settings,
            tab_log,
            tab_nexus,
        ]  # Add the new tab
        self.tabs.set_title(0, "Reduce")
        self.tabs.set_title(1, "Settings")
        self.tabs.set_title(2, "Log")
        self.tabs.set_title(3, "Nexus Explorer")  # Set the title for the new tab

        self.main = widgets.VBox(
            [
                self.proposal_number_box,
                self.tabs,
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

    def log_text(self, message):
        out = widgets.Output()
        with out:
            display(message)
        self.text_log.children = (out, *self.text_log.children)

    def log_progress(self, progress):
        self.progress_log.children = (progress,)

    def log_plot(self, plot):
        """Log a plot with a remove button and comment box."""
        # Create unique ID for this plot group - currently unused
        self.plot_counter += 1

        # Convert any existing interactive plots to static
        # NOTE: This conversion is not working as intended - plots remain interactive
        for child in self.plot_log.children:
            if isinstance(child, widgets.VBox):
                plot_output = child.children[1]  # Get the plot output widget
                with plot_output:
                    # Clear the output and redisplay as static
                    plot_output.clear_output()
                    if hasattr(plot_output, 'current_plot'):
                        display(plot_output.current_plot)

        # Create the plot output for the new plot
        plot_output = widgets.Output()
        plot_output.current_plot = plot

        # Create legend checkboxes container
        legend_container = widgets.HBox(
            layout=widgets.Layout(flex_wrap='wrap', align_items='center', padding='5px')
        )

        # Store original data for x4 transform
        original_data = {}

        # Store original methods before defining custom ones
        original_zoom = plot.view.canvas.zoom
        original_panzoom = plot.view.canvas.panzoom
        original_draw = plot.view.canvas.draw

        # Define all helper functions first
        def calculate_plot_limits(current_plot, is_transformed=False):
            """Calculate proper limits from all artists' data"""
            artists = current_plot.artists
            all_y_values = []
            all_y_errors = []
            for artist in artists.values():
                valid_mask = ~np.isnan(artist._data.values) & ~np.isinf(
                    artist._data.values
                )
                all_y_values.extend(artist._data.values[valid_mask])
                if artist._data.variances is not None:
                    error_mask = valid_mask & ~np.isinf(np.sqrt(artist._data.variances))
                    all_y_errors.extend(np.sqrt(artist._data.variances[error_mask]))

            if all_y_values:
                y_min = min(all_y_values)
                y_max = max(all_y_values)
                if all_y_errors:
                    y_min = min(y_min, min(all_y_values - np.array(all_y_errors)))
                    y_max = max(y_max, max(all_y_values + np.array(all_y_errors)))

                # Handle negative values for log scale
                if y_min <= 0:
                    positive_values = [y for y in all_y_values if y > 0]
                    if positive_values:
                        y_min = min(positive_values)
                        if all_y_errors:
                            positive_indices = [
                                i for i, y in enumerate(all_y_values) if y > 0
                            ]
                            error_values = [all_y_errors[i] for i in positive_indices]
                            if (
                                error_values
                            ):  # Only process if we have valid error values
                                y_min = min(
                                    y_min, min(positive_values - np.array(error_values))
                                )

                # Add padding (5% on log scale)
                log_range = np.log10(y_max) - np.log10(y_min)
                padding = 0.05 * log_range
                y_min = y_min * 10 ** (-padding)
                y_max = y_max * 10**padding

                # Use different minimum limits for transformed vs untransformed data
                if is_transformed:
                    # For y*x^4 data, use 3 orders of magnitude range
                    y_min = max(y_min, y_max * 1e-3)
                else:
                    # For original data, use 6 orders of magnitude range
                    y_min = max(y_min, y_max * 1e-6)

                return min(y_min, y_max), max(y_min, y_max)
            return None

        def fix_axis_orientation(current_plot, is_transformed=False):
            """Helper function to ensure correct axis orientation"""
            canvas = current_plot.view.canvas
            if hasattr(canvas, 'ax'):
                ax = canvas.ax

                # Calculate proper limits from data
                limits = calculate_plot_limits(current_plot, is_transformed)
                if limits is not None:
                    plot_y_min, plot_y_max = limits
                else:
                    # If no valid limits, use current ones
                    y_min, y_max = ax.get_ylim()
                    plot_y_min = min(y_min, y_max)
                    plot_y_max = max(y_min, y_max)

                # Force log scale and correct orientation
                ax.set_yscale('log')
                ax.set_ylim(plot_y_min, plot_y_max)
                if ax.yaxis.get_inverted():
                    ax.invert_yaxis()

                # Update canvas state to match
                canvas.yscale = 'log'
                canvas.yrange = (plot_y_min, plot_y_max)

                # Force a redraw to ensure changes take effect
                ax.figure.canvas.draw()

        def custom_reset_mode():
            """Custom reset mode that maintains correct orientation"""
            fix_axis_orientation(plot)

        def custom_zoom(event=None):
            """Custom zoom that maintains correct orientation"""
            if event is not None:
                original_zoom(event)
            else:
                original_zoom()
            fix_axis_orientation(plot)

        def custom_panzoom(event=None):
            """Custom panzoom that maintains correct orientation"""
            if event is not None:
                original_panzoom(event)
            else:
                original_panzoom()
            fix_axis_orientation(plot)

        def custom_draw(*args, **kwargs):
            """Custom draw that maintains correct orientation"""
            original_draw(*args, **kwargs)
            fix_axis_orientation(plot)

        def apply_x4_transform(change):
            """Apply or revert x4 transformation."""
            # Get the current plot
            current_plot = plot_output.current_plot
            # Get the figure's artists
            artists = current_plot.artists

            # Apply or revert transformation for each artist
            for name, artist in artists.items():
                if change['new']:  # Button is pressed - apply transformation
                    # Store original data if not already stored
                    if name not in original_data:
                        original_data[name] = artist._data.copy()

                    # Get x and y data from the artist
                    x_data = artist._coord.values
                    y_data = artist._data

                    # Make sure x_data matches y_data shape
                    if x_data.shape != y_data.values.shape:
                        x_centers = (x_data[1:] + x_data[:-1]) / 2
                        new_values = y_data.copy()
                        valid_mask = ~np.isnan(y_data.values)
                        transformed_values = y_data.values * x_centers**4
                        new_values.values = np.where(
                            valid_mask, transformed_values, y_data.values
                        )
                        if y_data.variances is not None:
                            new_values.variances = np.where(
                                valid_mask,
                                y_data.variances * x_centers**8,
                                y_data.variances,
                            )
                    else:
                        new_values = y_data.copy()
                        valid_mask = ~np.isnan(y_data.values)
                        transformed_values = y_data.values * x_data**4
                        new_values.values = np.where(
                            valid_mask, transformed_values, y_data.values
                        )
                        if y_data.variances is not None:
                            new_values.variances = np.where(
                                valid_mask,
                                y_data.variances * x_data**8,
                                y_data.variances,
                            )
                else:  # Button is unpressed - revert to original
                    new_values = original_data[name].copy()

                # Update the line with new values
                artist.update(new_values)

            # Update view limits to show all data
            fix_axis_orientation(current_plot, is_transformed=change['new'])

            # Update without autoscaling
            current_plot.update()

        # Now override the methods
        plot.view.canvas.reset_mode = custom_reset_mode
        plot.view.canvas.zoom = custom_zoom
        plot.view.canvas.panzoom = custom_panzoom
        plot.view.canvas.draw = custom_draw

        # Create and populate legend controls
        def create_legend_controls():
            legend_controls = []
            if hasattr(plot.view.canvas, 'ax'):
                ax = plot.view.canvas.ax
                handles, labels = ax.get_legend_handles_labels()

                # FIXME: Current limitation - Toggle does not affect error bars
                # This is because the error bars are separate artists in matplotlib
                # A fix would require finding and toggling associated error bar artists
                for handle, label in zip(handles, labels, strict=True):
                    checkbox = widgets.Checkbox(
                        value=True,
                        description=label,
                        indent=False,
                        layout=widgets.Layout(margin='0 10px 0 0'),
                    )

                    def make_toggle_callback(artist_handle):
                        def toggle_visibility(change):
                            artist_handle.set_visible(change['new'])
                            plot.view.canvas.draw()

                        return toggle_visibility

                    checkbox.observe(make_toggle_callback(handle), names='value')
                    legend_controls.append(checkbox)
            return legend_controls

        # Update legend container
        def update_legend():
            legend_container.children = create_legend_controls()

        # Display plot and create initial legend controls
        with plot_output:
            display(plot)
            update_legend()

        # Create buttons and comment box
        remove_button = widgets.Button(description='Remove plot')
        x4_button = widgets.ToggleButton(description='R*Q⁴')
        comment_box = widgets.Textarea(
            placeholder='Add comments about this plot here...',
            layout=widgets.Layout(width='75%', height='40px'),
        )

        # Connect the x4 transform button
        x4_button.observe(apply_x4_transform, names='value')

        # Create container for plot controls with legend
        controls_container = widgets.VBox(
            [  # Change to VBox for vertical stacking
                widgets.HBox(
                    [  # First row with buttons and comment
                        widgets.VBox(
                            [remove_button, x4_button]
                        ),  # Stack these vertically
                        comment_box,
                    ],
                    layout=widgets.Layout(width='100%'),
                ),
                widgets.VBox(
                    [  # Second row with legend controls
                        widgets.Label('Toggle Datasets:'),
                        legend_container,
                    ],
                    layout=widgets.Layout(margin='10px 0'),
                ),  # Add vertical margin
            ],
            layout=widgets.Layout(width='100%', align_items='flex-start'),
        )

        # Create vertical container for all elements
        plot_container = widgets.VBox([controls_container, plot_output])

        def remove_plot(_):
            # Find and remove this plot container
            self.plot_log.children = tuple(
                child for child in self.plot_log.children if child is not plot_container
            )

        remove_button.on_click(remove_plot)

        # Add the new plot container at the top
        self.plot_log.children = (plot_container, *self.plot_log.children)

    def create_hdf5_tree(self, filepath):
        """Create a tree representation of an HDF5 file structure."""

        def create_node(name, obj, path=''):
            full_path = f"{path}/{name}" if path else name
            if isinstance(obj, h5py.Dataset):
                # For datasets, show shape and dtype
                display_name = f"{name} ({obj.shape}, {obj.dtype})"
                node = Node(display_name, opened=False, icon='file')
                node.nexus_path = full_path  # Store path as custom attribute
                return node
            else:
                # For groups, create parent node and add children
                parent = Node(name, opened=False, icon='folder')
                parent.nexus_path = full_path  # Store path as custom attribute
                # Just iterate over the keys directly
                for child_name in obj.keys():
                    parent.add_node(create_node(child_name, obj[child_name], full_path))
                return parent

        try:
            with h5py.File(filepath, 'r') as f:
                root_node = create_node('', f)
                return Tree(nodes=[root_node])
        except Exception as e:
            # Use explicit conversion flag
            return Tree(nodes=[Node(f"Error loading file: {e!s}")])

    def update_nexus_view(self, *_):
        """Update the Nexus file viewer based on selected run."""
        selections = self.runs_table.selections
        if not selections:
            self.nexus_tree.nodes = [Node("Select a run to view its structure")]
            return

        # Get the first selected row
        row_idx = selections[0]['r1']
        run = self.runs_table.data.iloc[row_idx]['Run']
        filepath = self.get_filepath_from_run(run)

        # Create and display the tree for this file
        new_tree = self.create_hdf5_tree(filepath)
        self.nexus_tree.nodes = new_tree.nodes

    def display_nexus_content(self, path, h5file):
        """Display the content of a Nexus entry."""
        try:
            item = h5file[path] if path else h5file
            content = []

            # Show attributes if any
            if len(item.attrs) > 0:
                content.append("Attributes:")
                for key, value in item.attrs.items():
                    content.append(f"  {key}: {value}")

            # Show dataset content if it's a dataset
            if isinstance(item, h5py.Dataset):
                content.append("\nDataset content:")
                data = item[()]
                if data.size > 100:  # Truncate large datasets
                    content.append(f"  Shape: {data.shape}")
                    content.append("  First few values:")
                    content.append(f"  {data.flatten()[:100]}")
                    content.append("  ...")
                else:
                    content.append(f"  {data}")

            self.nexus_content.value = '\n'.join(content)
        except Exception as e:
            # Use explicit conversion flag
            self.nexus_content.value = f"Error reading content: {e!s}"

    def on_tree_select(self, event):
        """Handle tree node selection."""
        if not event['new']:  # No selection
            self.nexus_content.value = "Select a node to view its content"
            return

        selected_node = event['new'][0]

        # Get the path from the custom attribute
        path = getattr(selected_node, 'nexus_path', selected_node.name)

        # Get the currently selected run
        selections = self.runs_table.selections
        if not selections:
            return

        row_idx = selections[0]['r1']
        run = self.runs_table.data.iloc[row_idx]['Run']
        filepath = self.get_filepath_from_run(run)

        with h5py.File(filepath, 'r') as f:
            self.display_nexus_content(path, f)


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
        self._setdefault(df, "Comment", "")  # Add default empty comment
        df = self._ordercolumns(
            df, 'Run', 'Sample', 'Angle', 'Exclude', 'Comment'
        )  # Reorder columns
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

        def get_unique_names(df):
            # Create labels with Sample name and runs
            labels = [
                f"{params['Sample']} ({','.join(params['Runs'])})"
                for (_, params) in df.iterrows()
            ]
            duplicated_name_counter = {}
            unique = []
            for i, name in enumerate(labels):
                if name not in labels[:i]:
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
            figsize=(12, 6),  # Make figure wider initially
        )

        # Adjust the figure and legend
        if hasattr(all_runs_plot.view.canvas, 'ax'):
            ax = all_runs_plot.view.canvas.ax
            # Move legend outside
            ax.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
            # Adjust layout to prevent legend cutoff
            if hasattr(all_runs_plot.view.canvas, 'fig'):
                fig = all_runs_plot.view.canvas.fig
                fig.tight_layout()
                # Adjust subplot parameters to make room for legend
                fig.subplots_adjust(right=0.85)

        if matplotlib.get_backend().lower().startswith('qt'):
            all_runs_plot.show()
        else:
            self.log_plot(all_runs_plot)

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
        self.log_progress(progress)

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
