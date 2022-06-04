from __future__ import absolute_import
import sys
import os.path as op
from six.moves import range
sys.path.append('.')

from argparse import ArgumentParser, RawDescriptionHelpFormatter
from sfepy.base.base import nm, output
from sfepy.base.ioutils import remove_files
from sfepy.discrete.fem import Mesh, FEDomain
from sfepy.discrete.fem.meshio import output_mesh_formats, MeshIO
from sfepy.discrete.fem.mesh import fix_double_nodes
from sfepy.mesh.mesh_tools import triangulate

def convert_2d(filename_in = '1.mesh',filename_out = '1.h5',force_2d = True):

    mesh = Mesh.from_file(filename_in)

    if force_2d:
        data = list(mesh._get_io_data(cell_dim_only=2))
        data[0] = data[0][:, :2]
        mesh = Mesh.from_data(mesh.name, *data)

    output('writing %s...' % filename_out)
    mesh.write(filename_out, file_format=None, binary=False)
    output('...done')
