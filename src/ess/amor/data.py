# SPDX-License-Identifier: BSD-3-Clause
# Copyright (c) 2023 Scipp contributors (https://github.com/scipp)
import os

from ..reflectometry.types import Filename, ReferenceRun, SampleRun

_version = "1"


def _make_pooch():
    import pooch

    return pooch.create(
        path=pooch.os_cache("ess/amor"),
        env="ESS_AMOR_DATA_DIR",
        base_url="https://public.esss.dk/groups/scipp/ess/amor/{version}/",
        version=_version,
        registry={
            "reference.nxs": "md5:56d493c8051e1c5c86fb7a95f8ec643b",
            "sample.nxs": "md5:4e07ccc87b5c6549e190bc372c298e83",
            # Amor NeXus files
            # 608 to 611 are measurements on the same sample at different
            # sample rotations.
            # 612 and 613 are unknown to me.
            # 614 is a reference measurement on a super mirror.
            # The sample rotation values written in the files are wrong,
            # the real values are the following (all expressed in deg):
            #   608: 0.85,
            #   609: 2.25,
            #   610: 3.65,
            #   611: 5.05,
            #   612: 0.65,
            #   613: 0.65,
            #   614: 0.65,
            # The chopper phase offset value written in the files is wrong
            # the real value is -7.5 deg
            "amor2023n000608.hdf": "md5:e3a8b5b8495bb9ab2173848096df49d6",
            "amor2023n000609.hdf": "md5:d899d65f684cade2905f203f7d0fb326",
            "amor2023n000610.hdf": "md5:c9367d49079edcd17fa0b98e33326b05",
            "amor2023n000611.hdf": "md5:9da41177269faac0d936806393427837",
            "amor2023n000612.hdf": "md5:602f1bfcdbc1f618133c93a117d05f12",
            "amor2023n000613.hdf": "md5:ba0fbcbf0b45756a269eb3e943371ced",
            "amor2023n000614.hdf": "md5:18e8a755d6fd671758fe726de058e707",
            # Reflectivity curves obtained by applying Jochens Amor
            # software @ https://github.com/jochenstahn/amor.git
            # (repo commit hash e05fc9e1e124965919647f1856dbb9eb04221f1e).
            # to the above files.
            "608.Rqz.ort": "md5:60d8467796800f19c3aa1b6af5ad7b3d",
            "609.Rqz.ort": "md5:99af745e025423af64bc5a124a011826",
            "610.Rqz.ort": "md5:fa703dc9c5eed49d35e7a3d76a5746b9",
            "611.Rqz.ort": "md5:d6ea74d0525308a6a42938c106b62919",
            "612.Rqz.ort": "md5:f0b0eb269614645b4c6bbe24e29a6e10",
            "613.Rqz.ort": "md5:fff003552908c6a4a94b6213cfa08ac3",
        },
    )


def _make_pooch_protected():
    import pooch

    return pooch.create(
        path=pooch.os_cache("ess/protected/amor"),
        env="ESS_AMOR_DATA_DIR_PROTECTED",
        base_url="https://public.esss.dk/groups/scipp/protected/amor/{version}/",
        version=_version,
        registry={
            "test.md": "md5:d41d8cd98f00b204e9800998ecf8427e",
        },
    )


_pooch = _make_pooch()
_pooch_protected = _make_pooch_protected()


def amor_get_test():
    import base64

    import pooch

    username = os.environ.get('ESS_PROTECTED_FILESTORE_USERNAME', '')
    password = os.environ.get('ESS_PROTECTED_FILESTORE_PASSWORD', '')
    if username == '':
        raise RuntimeError('No credentials were found')

    return _pooch_protected.fetch(
        'test.md',
        downloader=pooch.HTTPDownloader(
            headers={
                'Authorization': 'Basic '
                + str(base64.b64encode(f'{username}:{password}'.encode()), 'utf-8')
            }
        ),
    )


def amor_old_sample_run() -> Filename[SampleRun]:
    return Filename[SampleRun](_pooch.fetch("sample.nxs"))


def amor_old_reference_run() -> Filename[ReferenceRun]:
    return Filename[ReferenceRun](_pooch.fetch("reference.nxs"))


def amor_reference_run() -> Filename[ReferenceRun]:
    return Filename[ReferenceRun](_pooch.fetch("amor2023n000614.hdf"))


def amor_sample_run(number: int | str) -> Filename[SampleRun]:
    return Filename[SampleRun](_pooch.fetch(f"amor2023n{int(number):06d}.hdf"))


def amor_psi_software_result(number: int | str) -> Filename[SampleRun]:
    return Filename[SampleRun](_pooch.fetch(f"{int(number):03d}.Rqz.ort"))


__all__ = [
    "amor_reference_run",
    "amor_sample_run",
    "amor_psi_software_result",
]
