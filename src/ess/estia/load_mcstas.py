# flake8: noqa
import numpy as np
import scipp as sc


def read_metadata(lines):
    '''Reads metadata of a McStas simulation from a .sim file'''
    data = {}
    section = None
    for line in lines:
        if line.startswith('begin'):
            _, _, name = line.partition(' ')
            section = {}
        elif line.startswith('end'):
            data.setdefault(name.strip(), []).append(section)
            section = None
        else:
            if section is not None:
                key, _, value = line.partition(': ')
                section[key.strip()] = value.strip()
    return data


def read_detector(lines):
    '''
    Reads detector event data from a McStas simulation.
    Assumes the first column is probabilities.
    '''
    meta = {}
    data = []
    for line in lines:
        if line.startswith('#'):
            key, _, value = line[2:].partition(': ')
            if '=' in value:
                key, _, value = value.partition('=')
            meta[key] = value
        else:
            data.append(list(map(float, line.strip().split(' '))))

    data = np.array(data)
    labels = meta['ylabel'].strip().split(' ')
    da = sc.DataArray(
        sc.array(dims=['events'], values=data[:, 0], variances=data[:, 0] ** 2),
        coords={
            label: sc.array(dims=['events'], values=values)
            for values, label in zip(data[:, 1:].T, labels[1:])
        },
    )
    for k, v in meta.items():
        da.coords[k] = sc.scalar(v)
    return da


def convert_mcstas_events_to_standard_events(da):
    da.coords["sample_rotation"] = sc.scalar(
        float(da.coords['omegaa'].value), unit='deg'
    ).to(unit='rad')
    da.coords["detector_rotation"] = 2 * da.coords["sample_rotation"]
    # da.coords["chopper_phase"] = chopper_phase
    da.coords["chopper_frequency"] = sc.scalar(14, unit='Hz') / sc.scalar(
        float(da.coords['enable_chopper'])
    )
    da.coords["L1"] = sc.scalar(39, unit='m')
    da.coords["L2"] = sc.scalar(4, unit='m')
    da.coords["angle_to_center_of_beam"] = angle_to_center_of_beam.to(unit='rad')
