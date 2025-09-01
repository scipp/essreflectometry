# ruff: noqa
from copy import deepcopy
from typing import TypeVar

import numpy as np
import scipp as sc
from scippnexus import NXdetector

import ess.estia as estia
import ess.reflectometry as reflectometry
from ess.amor.types import *
from ess.reduce.streaming import Accumulator, EternalAccumulator, StreamProcessor
from ess.reflectometry.types import *


def add_dummmy_value_to_all_transformations():
    import h5py

    with h5py.File('estia.hdf', 'r+') as f:

        def visitor(name, g):
            if (
                isinstance(g, h5py.Group)
                and g.attrs['NX_class'] == 'NXlog'
                and 'transformation_type' in g.attrs
            ):
                print(g)

                if len(g['time']) == 0:
                    print('Fixed')

                    attrs = dict(f[name + '/time'].attrs)
                    del f[name + '/time']
                    f.create_dataset(name + '/time', data=np.zeros(1))
                    for k, v in attrs.items():
                        f[name + '/time'].attrs[k] = v

                    attrs = dict(f[name + '/value'].attrs)
                    del f[name + '/value']
                    f.create_dataset(name + '/value', data=np.zeros(1))
                    for k, v in attrs.items():
                        f[name + '/value'].attrs[k] = v

        f.visititems(visitor)

    with h5py.File('estia.hdf', 'r') as f:

        def visitor(name, g):
            if (
                isinstance(g, h5py.Group)
                and g.attrs['NX_class'] == 'NXlog'
                and 'transformation_type' in g.attrs
            ):
                print(name, len(g['time']))

        f.visititems(visitor)


# add_dummmy_value_to_all_transformations()


def estia_stream_workflow() -> StreamProcessor:
    workflow = estia.EstiaWorkflow()
    workflow[SampleSize[SampleRun]] = sc.scalar(10.0, unit='mm')
    workflow[SampleSize[ReferenceRun]] = sc.scalar(10.0, unit='mm')

    workflow[ChopperPhase[ReferenceRun]] = sc.scalar(-7.5, unit='deg')
    workflow[ChopperPhase[SampleRun]] = sc.scalar(-7.5, unit='deg')

    workflow[QBins] = sc.geomspace(
        dim='Q', start=0.005, stop=0.3, num=391, unit='1/angstrom'
    )
    workflow[WavelengthBins] = sc.geomspace(
        'wavelength', 2.8, 12.5, 2001, unit='angstrom'
    )

    workflow[YIndexLimits] = sc.scalar(11), sc.scalar(41)
    workflow[ZIndexLimits] = sc.scalar(80), sc.scalar(370)
    workflow[BeamDivergenceLimits] = (
        sc.scalar(-0.75, unit='deg'),
        sc.scalar(0.75, unit='deg'),
    )

    workflow[SampleRotationOffset[SampleRun]] = sc.scalar(0.05, unit='deg')

    workflow[SampleRotationOffset[ReferenceRun]] = sc.scalar(0.05, unit='deg')

    workflow[Filename[SampleRun]] = './estia.hdf'
    workflow[Filename[ReferenceRun]] = './estia.hdf'

    workflow[ProtonCurrent[SampleRun]] = sc.DataArray(
        sc.ones(dims=['time'], shape=(1,)),
        coords={'time': sc.zeros(dims=['time'], shape=(1,), unit='s')},
    )
    workflow[ProtonCurrent[ReferenceRun]] = sc.DataArray(
        sc.ones(dims=['time'], shape=(1,)),
        coords={'time': sc.zeros(dims=['time'], shape=(1,), unit='s')},
    )
    workflow[DetectorSpatialResolution[SampleRun]] = 0.0025 * sc.units.m
    workflow[DetectorSpatialResolution[ReferenceRun]] = 0.0025 * sc.units.m
    import ess.reflectometry.supermirror as sm

    workflow[sm.CriticalEdge] = 0.022 * sc.Unit("1/angstrom")
    workflow[sm.Alpha] = sc.scalar(0.25 / 0.088, unit=sc.units.angstrom)
    workflow[sm.MValue] = 5

    workflow[SampleRotation[SampleRun]] = sc.scalar(1.0, unit='deg')
    workflow[DetectorRotation[SampleRun]] = sc.scalar(2.0, unit='deg')
    workflow[SampleRotation[ReferenceRun]] = sc.scalar(1.0, unit='deg')
    workflow[DetectorRotation[ReferenceRun]] = sc.scalar(2.0, unit='deg')

    workflow[DetectorData[ReferenceRun]] = workflow.compute(DetectorData[SampleRun])

    wf = workflow.copy()
    return StreamProcessor(
        wf,
        dynamic_keys=(NeXusComponent[NXdetector, SampleRun],),
        target_keys=(ReflectivityOverQ,),
        accumulators={
            Sample: EternalAccumulator,
            Reference: LatestAccumulator,
        },
    )


T = TypeVar('T')


class LatestAccumulator(Accumulator[T]):
    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self._value: T | None = None

    @property
    def is_empty(self) -> bool:
        return self._value is None

    def _get_value(self) -> T:
        return deepcopy(self._value)

    def _do_push(self, value: T) -> None:
        self._value = deepcopy(value)

    def clear(self) -> None:
        """Clear the accumulated value."""
        self._value = None


fake = sc.DataGroup(
    {
        'data': sc.DataArray(
            sc.ones(dims=['events'], shape=[10000]),
            coords={
                'event_id': sc.arange('events', 1, 10001),
                'event_time_offset': sc.arange('events', 1, 10001, unit='us'),
            },
        )
    }
)

s = estia_stream_workflow()
s.accumulate(
    {
        NeXusComponent[NXdetector, SampleRun]: fake,
    }
)
s.finalize()
