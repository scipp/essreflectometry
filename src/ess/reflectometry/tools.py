# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2025 Scipp contributors (https://github.com/scipp)
from __future__ import annotations

import uuid
from collections.abc import Callable, Mapping, Sequence
from itertools import chain
from typing import Any, Hashable

import numpy as np
import pandas as pd
import sciline as sl
import scipp as sc
import scipy.optimize as opt
from sciline.pipeline import _is_multiple_keys

from ess.reflectometry.types import (
    QBins,
    ReferenceRun,
    ReflectivityOverQ,
    SampleRotation,
    SampleRun,
    ScalingFactorForOverlap,
    UnscaledReducibleData,
)

_STD_TO_FWHM = sc.scalar(2.0) * sc.sqrt(sc.scalar(2.0) * sc.log(sc.scalar(2.0)))


def fwhm_to_std(fwhm: sc.Variable) -> sc.Variable:
    """
    Convert from full-width half maximum to standard deviation.

    Parameters
    ----------
    fwhm:
        Full-width half maximum.

    Returns
    -------
    :
        Standard deviation.
    """
    # Enables the conversion from full width half
    # maximum to standard deviation
    return fwhm / _STD_TO_FWHM


def linlogspace(
    dim: str,
    edges: list | np.ndarray,
    scale: list | str,
    num: list | int,
    unit: str | None = None,
) -> sc.Variable:
    """
    Generate a 1d array of bin edges with a mixture of linear and/or logarithmic
    spacings.

    Examples:

    - Create linearly spaced edges (equivalent to `sc.linspace`):
        linlogspace(dim='x', edges=[0.008, 0.08], scale='linear', num=50, unit='m')
    - Create logarithmically spaced edges (equivalent to `sc.geomspace`):
        linlogspace(dim='x', edges=[0.008, 0.08], scale='log', num=50, unit='m')
    - Create edges with a linear and a logarithmic part:
        linlogspace(dim='x', edges=[1, 3, 8], scale=['linear', 'log'], num=[16, 20])

    Parameters
    ----------
    dim:
        The dimension of the output Variable.
    edges:
        The edges for the different parts of the mesh.
    scale:
        A string or list of strings specifying the scaling for the different
        parts of the mesh. Possible values for the scaling are `"linear"` and `"log"`.
        If a list is supplied, the length of the list must be one less than the length
        of the `edges` parameter.
    num:
        An integer or a list of integers specifying the number of points to use
        in each part of the mesh. If a list is supplied, the length of the list must be
        one less than the length of the `edges` parameter.
    unit:
        The unit of the output Variable.

    Returns
    -------
    :
        Lin-log spaced Q-bin edges.
    """
    if not isinstance(scale, list):
        scale = [scale]
    if not isinstance(num, list):
        num = [num]
    if len(scale) != len(edges) - 1:
        raise ValueError(
            "Sizes do not match. The length of edges should be one greater than scale."
        )

    funcs = {"linear": sc.linspace, "log": sc.geomspace}
    grids = []
    for i in range(len(edges) - 1):
        # Skip the leading edge in the piece when concatenating
        start = int(i > 0)
        mesh = funcs[scale[i]](
            dim=dim, start=edges[i], stop=edges[i + 1], num=num[i] + start, unit=unit
        )
        grids.append(mesh[dim, start:])

    return sc.concat(grids, dim)


class BatchWorkflow:
    """ """

    def __init__(
        self,
        workflow: sl.Pipeline,
        param_table: pd.DataFrame | None = None,
        # groupings: list[dict] | None = None,
    ):
        self.workflow = workflow.copy()
        if param_table is not None:
            self.workflow = self.workflow.map(param_table)

    #     # self._original_workflow = workflow.copy()
    #     # self.param_table = param_table.copy()
    #     # self._groupings = groupings.copy() if groupings is not None else []
    #     if self.param_table is not None:
    #         self.workflow = self._original_workflow.map(self.param_table)
    #     else:
    #         self._mapped_workflow = self._original_workflow.copy()
    #     # self._mapped_workflow = workflow.map(param_table)

    # # @classmethod
    # # def _from_mapped(cls, mapped_workflow) -> 'WorkflowCollection':
    # #     out = cls(sl.Pipeline(), {})
    # #     out._mapped_workflow = mapped_workflow.copy()
    # #     return out

    def __setitem__(self, key, value) -> None:
        # # TODO: the setitem is currently broken if the workflow has been grouped
        # if key in self.param_table:
        #     ind = list(self.param_table.keys()).index(key)
        #     self.param_table.iloc[:, ind] = value
        #     # self._mapped_workflow = self._original_workflow.map(self.param_table)
        # else:
        #     self._original_workflow[key] = value
        #     # self.param_table.insert(len(self.param_table.columns), key, value)
        #     # # self._original_workflow[key] = None
        #     # self._mapped_workflow = self._original_workflow.map(self.param_table)

        # Setting a value that is in the original parameter table of the mapped workflow
        # does not currently work (see https://github.com/scipp/sciline/issues/224)
        # So we have to go and modify the underlying parameter table directly
        if key in self.workflow._cbgraph._node_values._values:
            self.workflow._cbgraph._node_values._values[key]._series.update(value)
        else:
            if sl.is_mapped_node(self.workflow, key):
                targets = sl.get_mapped_node_names(self.workflow, key)
                for k, v in value.items():
                    self.workflow[targets[k]] = v
            else:
                self.workflow[key] = value

    def __getitem__(self, key) -> BatchWorkflow:
        return BatchWorkflow(
            workflow=self.workflow[key],
            # param_table=self.param_table,
            # groupings=self._groupings,
        )

    def compute(
        self,
        keys: type | Sequence[type],
        index_names: Sequence[Hashable] | None = None,
        **kwargs,
    ) -> Mapping[str, Any]:
        # from sciline.pipeline import _is_multiple_keys

        # mapped_workflow = self._build_workflow()

        out = {}
        if not _is_multiple_keys(keys):
            keys = [keys]
        for key in keys:
            out[key] = {}
            if sl.is_mapped_node(self.workflow, key):
                targets = sl.get_mapped_node_names(
                    self.workflow, key, index_names=index_names
                )
                results = self.workflow.compute(targets, **kwargs)
                for node, v in results.items():
                    out[key][node.index.values[0]] = v
            else:
                out[key] = self.workflow.compute(key, **kwargs)
        return next(iter(out.values())) if len(out) == 1 else out

    # TODO: implement get()

    # TODO: adding groupby lead to too many headaches: having to keep the original
    # workflow around so we could map over unique groups to re-attach the bottom part
    # of the graph after grouping, having to then keep track of successive groupings
    # because the workflow was being built from scratch each time, etc.

    # # def groupby(
    # #     self, *, at: str | type, key: str | type, reduce: Callable
    # # ) -> 'WorkflowCollection':
    # #     return WorkflowCollection(
    # #         workflow=self._original_workflow,
    # #         param_table=self.param_table,
    # #         groupings=[*self._groupings, {'at': at, 'key': key, 'reduce': reduce}],
    # #     )

    # def _apply_groupings(self) -> sl.Pipeline:
    #     # original_workflow = self._original_workflow.copy()

    #     # # Split the param table into the keys that are used for grouping and the rest
    #     # grouping_table = self.param_table[
    #     #     list({g['key'] for g in self._groupings})
    #     # ].copy()
    #     # mapping_table = self.param_table.drop(columns=grouping_table.columns)

    #     # TODO: should this only map the grouping table?
    #     # grouped = self._original_workflow.map(self.param_table)
    #     grouped = self._original_workflow.map(self.param_table)
    #     node_names = []
    #     for grouping in self._groupings:
    #         at = grouping['at']
    #         key = grouping['key']
    #         reduce = grouping['reduce']
    #         node_names.append(uuid.uuid4().hex)
    #         grouped = (
    #             # grouped.map(self.param_table)  # [at]
    #             grouped.groupby(key).reduce(key=at, func=reduce, name=node_names[-1])
    #         )

    #     # new = self._original_workflow.map(mapping_table)
    #     new = self._original_workflow.copy()
    #     for name, grouping in zip(node_names, self._groupings, strict=True):
    #         # new = self._original_workflow.copy()
    #         at = grouping['at']
    #         key = grouping['key']
    #         reduce = grouping['reduce']
    #         new[at] = None
    #         unique_groups = sl.get_mapped_node_names(grouped, name).index

    #         # Note that we convert the key to a string, because pandas DataFrame does
    #         # not support using types as index name.
    #         str_key = str(key)
    #         new = new.map(
    #             pd.DataFrame(
    #                 {at: [None] * len(unique_groups), str_key: unique_groups}
    #             ).set_index(str_key)
    #         )

    #     for name, grouping in zip(node_names, self._groupings, strict=True):
    #         # at = grouping['at']
    #         new[grouping['at']] = grouped[name]

    #     # # Build param table from the entries in the mapping table that are not yet
    #     # # mapped
    #     # missing_keys = [
    #     #     k for k in mapping_table.columns if not sl.is_mapped_node(new, k)
    #     # ]
    #     # print(f'Mapping missing keys: {missing_keys}')
    #     # if missing_keys:
    #     #     print("Mapping table", mapping_table[missing_keys])
    #     #     new = new.map(mapping_table[missing_keys])
    #     return new
    #     # return WorkflowCollection(workflow=mapped, param_table=None)

    # def _build_workflow(self) -> sl.Pipeline:
    #     if not self._groupings:
    #         return self._original_workflow.map(self.param_table)
    #     else:
    #         return self._apply_groupings()

    def visualize(
        self, targets, index_names: Sequence[Hashable] | None = None, **kwargs
    ):
        # mapped_workflow = self._build_workflow()
        if not _is_multiple_keys(targets):
            targets = [targets]
        types = []
        for t in targets:
            if sl.is_mapped_node(self.workflow, t):
                types.extend(
                    sl.get_mapped_node_names(self.workflow, t, index_names=index_names)
                )
            else:
                types.append(t)
        key = types[0] if len(types) == 1 else types
        # if sl.is_mapped_node(self._mapped_workflow, key):
        # targets = sl.get_mapped_node_names(self._mapped_workflow, targets)
        return self.workflow.visualize(key, **kwargs)

    def copy(self) -> BatchWorkflow:
        return BatchWorkflow(
            workflow=self.workflow,
            # param_table=self.param_table,
            # groupings=None if not self._groupings else self._groupings,
        )


def _sort_by(a, by):
    return [x for x, _ in sorted(zip(a, by, strict=True), key=lambda x: x[1])]


def _find_interval_overlaps(intervals):
    '''Returns the intervals where at least
    two or more of the provided intervals
    are overlapping.'''
    edges = list(chain.from_iterable(intervals))
    is_start_edge = list(chain.from_iterable((True, False) for _ in intervals))
    edges_sorted = sorted(edges)
    is_start_edge_sorted = _sort_by(is_start_edge, edges)

    number_overlapping = 0
    overlap_intervals = []
    for x, is_start in zip(edges_sorted, is_start_edge_sorted, strict=True):
        if number_overlapping == 1 and is_start:
            start = x
        if number_overlapping == 2 and not is_start:
            overlap_intervals.append((start, x))
        if is_start:
            number_overlapping += 1
        else:
            number_overlapping -= 1
    return overlap_intervals


def _searchsorted(a, v):
    for i, e in enumerate(a):
        if e > v:
            return i
    return len(a)


def _create_qgrid_where_overlapping(qgrids):
    '''Given a number of Q-grids, construct a new grid
    covering the regions where (any two of the) provided grids overlap.'''
    pieces = []
    for start, end in _find_interval_overlaps([(q.min(), q.max()) for q in qgrids]):
        interval_sliced_from_qgrids = [
            q[max(_searchsorted(q, start) - 1, 0) : _searchsorted(q, end) + 1]
            for q in qgrids
        ]
        densest_grid_in_interval = max(interval_sliced_from_qgrids, key=len)
        pieces.append(densest_grid_in_interval)
    return sc.concat(pieces, dim='Q')


def _same_dtype(arrays):
    return [arr.to(dtype='float64') for arr in arrays]


def _interpolate_on_qgrid(curves, grid):
    return sc.concat(
        _same_dtype([sc.lookup(c, grid.dim)[sc.midpoints(grid)] for c in curves]),
        dim='curves',
    )


def scale_reflectivity_curves_to_overlap(
    workflow: BatchWorkflow | sl.Pipeline,
    critical_edge_interval: tuple[sc.Variable, sc.Variable] | None = None,
    cache_intermediate_results: bool = True,
) -> tuple[list[sc.DataArray], list[sc.Variable]]:
    '''
    Set the ``ScalingFactorForOverlap`` parameter on the provided workflows
    in a way that would makes the 1D reflectivity curves overlap.
    One can supply either a collection of workflows or a single workflow.

    If :code:`critical_edge_interval` is not provided, all workflows are scaled except
    the data with the lowest Q-range, which is considered to be the reference curve.
    The scaling factors are determined by a maximum likelihood estimate
    (assuming the errors are normal distributed).

    If :code:`critical_edge_interval` is provided then all data are scaled.

    All reflectivity curves must be have the same unit for data and the Q-coordinate.

    Parameters
    ---------
    workflows:
        The workflow or collection of workflows that can compute ``ReflectivityOverQ``.
    critical_edge_interval:
        A tuple denoting an interval that is known to belong
        to the critical edge, i.e. where the reflectivity is
        known to be 1.
    cache_intermediate_results:
        If ``True`` the intermediate results ``UnscaledReducibleData`` will be cached
        (this is the base for all types that are downstream of the scaling factor).

    Returns
    ---------
    :
        A list of scaled reflectivity curves and a list of the scaling factors.
    '''
    not_batch = isinstance(workflow, sl.Pipeline)

    batch = workflow.copy()
    if cache_intermediate_results:
        try:
            batch[UnscaledReducibleData[SampleRun]] = batch.compute(
                UnscaledReducibleData[SampleRun]
            )
        except sl.UnsatisfiedRequirement:
            pass
        try:
            batch[UnscaledReducibleData[ReferenceRun]] = batch.compute(
                UnscaledReducibleData[ReferenceRun]
            )
        except sl.UnsatisfiedRequirement:
            pass

    reflectivities = batch.compute(ReflectivityOverQ)
    if not_batch:
        reflectivities = {"": reflectivities}

    # First sort the dict of reflectivities by the Q min value
    curves = {
        k: v.hist() if v.bins is not None else v
        for k, v in sorted(
            reflectivities.items(), key=lambda item: item[1].coords['Q'].min().value
        )
    }

    critical_edge_key = uuid.uuid4().hex
    if critical_edge_interval is not None:
        q = batch.compute(QBins)
        if hasattr(q, "items"):
            # If QBins is a mapping, find the one with the lowest Q start
            # Note the conversion to a dict, because if pandas is used for the mapping,
            # it will return a Series, whose `.values` attribute is not callable.
            q = min(dict(q).values(), key=lambda q_: q_.min())

        # TODO: This is slightly different from before: it extracts the bins from the
        # QBins variable that cover the critical edge interval. This means that the
        # resulting curve will not necessarily begin and end exactly at the values
        # specified, but rather at the closest bin edges.
        edge = sc.DataArray(
            data=sc.ones(sizes={q.dim: q.sizes[q.dim] - 1}, with_variances=True),
            coords={q.dim: q},
        )[q.dim, critical_edge_interval[0] : critical_edge_interval[1]]
        # Now place the critical edge at the beginning
        curves = {critical_edge_key: edge} | curves

    if len({c.data.unit for c in curves.values()}) != 1:
        raise ValueError('The reflectivity curves must have the same unit')
    if len({c.coords['Q'].unit for c in curves.values()}) != 1:
        raise ValueError('The Q-coordinates must have the same unit for each curve')

    qgrid = _create_qgrid_where_overlapping([c.coords['Q'] for c in curves.values()])

    r = _interpolate_on_qgrid(map(sc.values, curves.values()), qgrid).values
    v = _interpolate_on_qgrid(map(sc.variances, curves.values()), qgrid).values

    def cost(scaling_factors):
        scaling_factors = np.concatenate([[1.0], scaling_factors])[:, None]
        r_scaled = scaling_factors * r
        v_scaled = scaling_factors**2 * v
        v_scaled[v_scaled == 0] = np.nan
        inv_v_scaled = 1 / v_scaled
        r_avg = np.nansum(r_scaled * inv_v_scaled, axis=0) / np.nansum(
            inv_v_scaled, axis=0
        )
        return np.nansum((r_scaled - r_avg) ** 2 * inv_v_scaled)

    sol = opt.minimize(cost, [1.0] * (len(curves) - 1))
    scaling_factors = (1.0, *map(float, sol.x))

    results = {
        k: v
        for k, v in zip(curves.keys(), scaling_factors, strict=True)
        if k != critical_edge_key
    }
    if not_batch:
        results = results[""]
    batch[ScalingFactorForOverlap[SampleRun]] = results
    # return results

    return batch


def combine_curves(
    curves: Sequence[sc.DataArray],
    q_bin_edges: sc.Variable | None = None,
) -> sc.DataArray:
    '''Combines the given curves by interpolating them
    on a 1d grid defined by :code:`q_bin_edges` and averaging
    over the provided reflectivity curves.

    The averaging is done using a weighted mean where the weights
    are proportional to the variances.

    Unless the curves are already scaled correctly they might
    need to be scaled using :func:`scale_reflectivity_curves_to_overlap`
    before calling this function.

    All curves must be have the same unit for data and the Q-coordinate.

    Parameters
    ----------
    curves:
        the reflectivity curves that should be combined
    q_bin_edges:
        the Q bin edges of the resulting combined reflectivity curve

    Returns
    ---------
    :
        A data array representing the combined reflectivity curve
    '''
    if len({c.data.unit for c in curves}) != 1:
        raise ValueError('The reflectivity curves must have the same unit')
    if len({c.coords['Q'].unit for c in curves}) != 1:
        raise ValueError('The Q-coordinates must have the same unit for each curve')

    r = _interpolate_on_qgrid(map(sc.values, curves), q_bin_edges).values
    v = _interpolate_on_qgrid(map(sc.variances, curves), q_bin_edges).values

    v[v == 0] = np.nan
    inv_v = 1.0 / v
    r_avg = np.nansum(r * inv_v, axis=0) / np.nansum(inv_v, axis=0)
    v_avg = 1 / np.nansum(inv_v, axis=0)
    return sc.DataArray(
        data=sc.array(
            dims='Q',
            values=r_avg,
            variances=v_avg,
            unit=next(iter(curves)).data.unit,
        ),
        coords={'Q': q_bin_edges},
    )


class NoParameter:
    def __repr__(self):
        return "<NOPARAM>"


NO_PARAMETER = NoParameter()


def parameter_table(params: Mapping[Any, Mapping[type, Any]]) -> pd.DataFrame:
    """
    Create a parameter table from the provided params.

    Example:

    ```
    runs = {
        '608': {
            SampleRotationOffset[SampleRun]: sc.scalar(0.05, unit='deg'),
            Filename[SampleRun]: amor.data.amor_run(608),
        },
        '609': {
            SampleRotationOffset[SampleRun]: sc.scalar(0.06, unit='deg'),
            Filename[SampleRun]: amor.data.amor_run(609),
        },
        '610': {
            SampleRotationOffset[SampleRun]: sc.scalar(0.05, unit='deg'),
            Filename[SampleRun]: amor.data.amor_run(610),
        },
        '611': {
            SampleRotationOffset[SampleRun]: sc.scalar(0.07, unit='deg'),
            Filename[SampleRun]: amor.data.amor_run(611),
        },
    }

    param_table = parameter_table(runs)
    ```

    Parameters
    ----------
    params:
        The sciline parameters to be used for each run.
        Should be a mapping where the keys are the names of the runs
        and the values are mappings of type to value pairs.
    """
    all_types = {t for v in params.values() for t in v.keys()}
    data = {t: [] for t in all_types}
    for param in params.values():
        for t in all_types:
            if t in param:
                data[t].append(param[t])
            else:
                data[t].append(NO_PARAMETER)

    return pd.DataFrame(data, index=params.keys()).rename_axis(index='run_id')


def merge(*das):
    da = sc.reduce(das).bins.concat()
    missing_coords = set(das[0].coords) - set(da.coords)
    return da.assign_coords({coord: das[0].coords[coord] for coord in missing_coords})


def batch_processor(
    workflow: sl.Pipeline,
    param_table: pd.DataFrame,
    group_by_sample_rotation: bool = False,
) -> BatchWorkflow:
    """
    Maps the provided workflow over the provided params.

    Example:

    ```
    from ess.reflectometry import amor, tools

    workflow = amor.AmorWorkflow()

    runs = {
        '608': {
            SampleRotationOffset[SampleRun]: sc.scalar(0.05, unit='deg'),
            Filename[SampleRun]: amor.data.amor_run(608),
        },
        '609': {
            SampleRotationOffset[SampleRun]: sc.scalar(0.06, unit='deg'),
            Filename[SampleRun]: amor.data.amor_run(609),
        },
        '610': {
            SampleRotationOffset[SampleRun]: sc.scalar(0.05, unit='deg'),
            Filename[SampleRun]: amor.data.amor_run(610),
        },
        '611': {
            SampleRotationOffset[SampleRun]: sc.scalar(0.07, unit='deg'),
            Filename[SampleRun]: amor.data.amor_run(611),
        },
    }

    batch = tools.batch_processor(workflow, parameter_table(runs))

    results = batch.compute(ReflectivityOverQ)
    ```

    Parameters
    ----------
    workflow:
        The sciline workflow used to compute the targets for each of the runs.
    params:
        The sciline parameters to be used for each run.
        Should be a mapping where the keys are the names of the runs
        and the values are mappings of type to value pairs.
    """
    # import pandas as pd

    # all_types = {t for v in params.values() for t in v.keys()}
    # data = {t: [] for t in all_types}
    # for param in params.values():
    #     for t in all_types:
    #         if t in param:
    #             data[t].append(param[t])
    #         else:
    #             # Set the default value
    #             data[t].append(workflow.compute(t))

    # param_table = parameter_table(params)

    param_table = param_table.copy()

    # If ScaleingFactorForOverlap is not in the param table, add it
    if ScalingFactorForOverlap[SampleRun] not in param_table:
        param_table = param_table.join(
            pd.DataFrame(
                {ScalingFactorForOverlap[SampleRun]: [NO_PARAMETER] * len(param_table)},
                index=param_table.index,
            )
        )

    # Replace NO_PARAMETER entries by the default value
    for icol, tp in enumerate(param_table.columns):
        column = param_table.iloc[:, icol]
        if NO_PARAMETER in column.values:
            default = workflow.compute(tp)
            column[column.values == NO_PARAMETER] = default

    # TODO: If non-source nodes are in the param table, we need to remove them from the
    # mapping and override them in the workflow afterwards

    wf = BatchWorkflow(workflow, param_table)

    if not group_by_sample_rotation:
        return wf

    # key = SampleRotation[SampleRun]

    # Make a unique string key for the sample rotation. We don't use the
    # SampleRotation from the graph because it is not a source node

    at = UnscaledReducibleData[SampleRun]
    key = SampleRotation[SampleRun]  # uuid.uuid4().hex
    # Make a string key for the sample rotation. We don't use the
    # SampleRotation from the graph because it is not a source node
    str_key = str(key)

    if key in param_table:
        # If the sample rotation is already in the param table, we need to drop it
        # from the param table because it is not a source node. Instead we copy the
        # values from the param table to a new column with a unique name
        pass

    else:
        # grouped = self._original_workflow.map(self.param_table)
        # at = UnscaledReducibleData[SampleRun]
        rotations = wf.compute(SampleRotation[SampleRun])
        # scalings = wf.compute(ScalingFactorForOverlap[SampleRun])

        # Groupby does not like Scipp variables, so we convert to floats for grouping
        rotation_unit = next(iter(rotations.values())).unit
        # TODO: make sure that the order of rotations is the same as in the param table
        param_table = param_table.join(
            pd.DataFrame(
                {key: [r.to(unit=rotation_unit).value for r in rotations.values()]},
                index=param_table.index,
            )
        )

    # # Map again to add the sample rotation column
    # wf = BatchWorkflow(workflow, param_table)

    grouping_node_name = uuid.uuid4().hex

    grouped = (
        workflow.map(param_table)
        .groupby(key)
        .reduce(key=at, func=merge, name=grouping_node_name)
    )

    new = workflow.copy()
    new[at] = None

    unique_groups = sl.get_mapped_node_names(grouped, grouping_node_name).index

    # Note that we convert the key to a string, because pandas DataFrame does
    # not support using types as index name.
    str_key = str(key)
    scale = workflow.compute(ScalingFactorForOverlap[SampleRun])
    new = new.map(
        pd.DataFrame(
            {
                at: [None] * len(unique_groups),
                str_key: unique_groups,
                ScalingFactorForOverlap[SampleRun]: [scale] * len(unique_groups),
            }
        ).set_index(str_key)
    )
    new[at] = grouped[grouping_node_name]

    # new = self._original_workflow.copy()
    # for name, grouping in zip(node_names, self._groupings, strict=True):
    #     # new = self._original_workflow.copy()
    #     at = grouping['at']
    #     key = grouping['key']
    #     reduce = grouping['reduce']
    #     new[at] = None
    #     unique_groups = sl.get_mapped_node_names(grouped, name).index

    #     # Note that we convert the key to a string, because pandas DataFrame does
    #     # not support using types as index name.
    #     str_key = str(key)
    #     new = new.map(
    #         pd.DataFrame(
    #             {at: [None] * len(unique_groups), str_key: unique_groups}
    #         ).set_index(str_key)
    #     )

    # for name, grouping in zip(node_names, self._groupings, strict=True):
    #     # at = grouping['at']
    #     new[grouping['at']] = grouped[name]

    # # # Build param table from the entries in the mapping table that are not yet
    # # # mapped
    # # missing_keys = [
    # #     k for k in mapping_table.columns if not sl.is_mapped_node(new, k)
    # # ]
    # # print(f'Mapping missing keys: {missing_keys}')
    # # if missing_keys:
    # #     print("Mapping table", mapping_table[missing_keys])
    # #     new = new.map(mapping_table[missing_keys])
    # return new
    wf = BatchWorkflow(new)

    return wf
