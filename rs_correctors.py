#!/usr/bin/env python
"""
Compute homogenized elastic coefficients for a given reconstructed microstructure of porous media or a composite material.
"""
from __future__ import print_function
from __future__ import absolute_import
from argparse import ArgumentParser
import sys
import six
sys.path.append('.')

import numpy as nm
import os

from sfepy import data_dir
import sfepy.discrete.fem.periodic as per
from sfepy.homogenization.utils import define_box_regions
import sfepy.mechanics.matcoefs as sff

def get_bord(coors, domain=None):
    x, y = coors[:, 0], coors[:, 1]
    r = nm.maximum(abs(x),abs(y))
    flag = nm.where((r > 0.998))[0]
    return flag

def define_regions():
    """
    Define various subdomains for a given mesh file.
    """
    regions = {}
    dim = 2

    regions['Y'] = 'all'
    regions['Ym'] = 'cells of group 1'
    regions['Yc'] = 'r.Y -c r.Ym'
    regions['Surf'] = ('vertices by get_bord','facet')
    lbp = (-1,-1)
    rtp = (1, 1)
    regions.update(define_box_regions(2, lbp, rtp))

    return dim, regions

##
# Mesh file.
filename_mesh = 'Input/meshfile.h5'

##
# Define regions (subdomains, boundaries) - $Y$, $Y_i$, ...
# depending on a mesh used.
dim, regions = define_regions()

functions = {
    'match_x_plane' : (per.match_x_plane,),
    'match_y_plane' : (per.match_y_plane,),
    'match_z_plane' : (per.match_z_plane,),
    'match_x_line' : (per.match_x_line,),
    'match_y_line' : (per.match_y_line,),
    'get_bord' : (get_bord,),
}

##
# Define fields: 'displacement' in $Y$,
# 'pressure_m' in $Y_m$.
fields = {
    'displacement' : ('real', dim, 'Y', 1),
}

##
# Define corrector variables: unknown displaements: uc, test: vc
# displacement-like variables: Pi, Pi1, Pi2
variables = {
    'uc'       : ('unknown field',   'displacement', 0),
    'vc'       : ('test field',      'displacement', 'uc'),
    'Pi'       : ('parameter field', 'displacement', 'uc'),
    'Pi1'      : ('parameter field', 'displacement', None),
    'Pi2'      : ('parameter field', 'displacement', None),
}

##
# Periodic boundary conditions.
if dim == 3:
    epbcs = {
        'periodic_x' : (['Left', 'Right'], {'uc.all' : 'uc.all'},
                        'match_x_plane'),
        'periodic_y' : (['Near', 'Far'], {'uc.all' : 'uc.all'},
                        'match_y_plane'),
        'periodic_z' : (['Top', 'Bottom'], {'uc.all' : 'uc.all'},
                        'match_z_plane'),
    }
else:
    epbcs = {
        'periodic_x': (['Left', 'Right'], {'uc.all': 'uc.all'},
                       'match_y_line'),
        'periodic_y': (['Bottom', 'Top'], {'uc.all': 'uc.all'},
                       'match_x_line'),
    }

##
# Dirichlet boundary conditions.
ebcs = {
    'fixed_u' : ('Corners', {'uc.all' : 0.0}),
}

##
# Material defining constitutive parameters of the microproblem.
materials = {
    'm' : ({'D' : {
        'Ym': sff.stiffness_from_youngpoisson(dim, 100, 0.2, plane='strain'),
        'Yc': sff.stiffness_from_youngpoisson(dim, 10^(-8), 10^(-8), plane='strain')}
    },),
}

##
# Numerical quadratures for volume (i3 - order 3) integral terms.
integrals = {
    'i3' : 3,
}

##
# Homogenized coefficients to compute.
def set_elastic(variables, ir, ic, mode, pis, corrs_rs):
    mode2var = {'row' : 'Pi1', 'col' : 'Pi2'}

    val = pis.states[ir, ic]['uc'] + corrs_rs.states[ir, ic]['uc']

    variables[mode2var[mode]].set_data(val)

coefs = {
    'E' : {
        'requires' : ['pis', 'corrs_rs'],
        'expression' : 'dw_lin_elastic.i3.Y(m.D, Pi1, Pi2)',
        'set_variables' : set_elastic,
    },
}

all_periodic = ['periodic_%s' % ii for ii in ['x', 'y', 'z'][:dim] ]
requirements = {
    'pis' : {
        'variables' : ['uc'],
    },
    ##
    # Steady state correctors $\bar{\omega}^{rs}$.
    'corrs_rs' : {
        'requires' : ['pis'],
        'save_variables' : ['uc'],
        'ebcs' : ['fixed_u'],
        'epbcs' : all_periodic,
        'equations' : {'eq' : """dw_lin_elastic.i3.Y(m.D, vc, uc)
                             = - dw_lin_elastic.i3.Y(m.D, vc, Pi)"""},
        'set_variables' : [('Pi', 'pis', 'uc')],
        'save_name' : 'corrs_elastic',
        'is_linear' : True,
    },
}

##
# Solvers.
solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max' : 1,
        'eps_a' : 1e-8,
        'eps_r' : 1e-2,
    })
}

############################################
# Mini-application below, computing the homogenized elastic coefficients.
helps = {
    'no_pauses' : 'do not make pauses',
}

def main():
    import os
    from sfepy.base.base import spause, output
    from sfepy.base.conf import ProblemConf, get_standard_keywords
    from sfepy.discrete import Problem
    import sfepy.homogenization.coefs_base as cb

    # parser = ArgumentParser(description=__doc__)
    # parser.add_argument('--version', action='version', version='%(prog)s')
    # parser.add_argument('-n', '--no-pauses',
    #                     action="store_true", dest='no_pauses',
    #                     default=False, help=helps['no_pauses'])
    # options = parser.parse_args()

    # if options.no_pauses:
    #     def spause(*args):
    #         output(*args)

    nm.set_printoptions(precision=3)

#     spause(r""">>>
# First, this file will be read in place of an input
# (problem description) file.
# Press 'q' to quit the example, press any other key to continue...""")
    required, other = get_standard_keywords()
    required.remove('equations')
    # Use this file as the input file.
    conf = ProblemConf.from_file(__file__, required, other)
    print(list(conf.to_dict().keys()))
#     spause(r""">>>
# ...the read input as a dict (keys only for brevity).
# ['q'/other key to quit/continue...]""")

#     spause(r""">>>
# Now the input will be used to create a Problem instance.
# ['q'/other key to quit/continue...]""")
    problem = Problem.from_conf(conf, init_equations=False)
    # The homogenization mini-apps need the output_dir.
    output_dir = ''
    problem.output_dir = output_dir
    print(problem)
#     spause(r""">>>
# ...the Problem instance.
# ['q'/other key to quit/continue...]""")

#     spause(r""">>>
# The homogenized elastic coefficient $E_{ijkl}$ is expressed
# using $\Pi$ operators, computed now. In fact, those operators are permuted
# coordinates of the mesh nodes.
# ['q'/other key to quit/continue...]""")
    req = conf.requirements['pis']
    mini_app = cb.ShapeDimDim('pis', problem, req)
    mini_app.setup_output(save_formats=['vtk'],
                          file_per_var=False)
    pis = mini_app()
    print(pis)
#     spause(r""">>>
# ...the $\Pi$ operators.
# ['q'/other key to quit/continue...]""")

#     spause(r""">>>
# Next, $E_{ijkl}$ needs so called steady state correctors $\bar{\omega}^{rs}$,
# computed now.
# ['q'/other key to quit/continue...]""")
    req = conf.requirements['corrs_rs']

    save_name = req.get('save_name', '')
    name = os.path.join(output_dir, save_name)

    mini_app = cb.CorrDimDim('steady rs correctors', problem, req)
    mini_app.setup_output(save_formats=['vtk'],
                          file_per_var=False)
    corrs_rs = mini_app(data={'pis': pis})
    print(corrs_rs)
#     spause(r""">>>
# ...the $\bar{\omega}^{rs}$ correctors.
# The results are saved in: %s.%s
#
# Try to display them with:
#
#    python postproc.py %s.%s
#
# ['q'/other key to quit/continue...]""" % (2 * (name, problem.output_format)))
#
#     spause(r""">>>
# Then the volume of the domain is needed.
# ['q'/other key to quit/continue...]""")
    volume = problem.evaluate('ev_volume.i3.Y(uc)')
    print(volume)

#     spause(r""">>>
# ...the volume.
# ['q'/other key to quit/continue...]""")
#
#     spause(r""">>>
# Finally, $E_{ijkl}$ can be computed.
# ['q'/other key to quit/continue...]""")
    mini_app = cb.CoefSymSym('homogenized elastic tensor',
                             problem, conf.coefs['E'])
    c_e = mini_app(volume, data={'pis': pis, 'corrs_rs' : corrs_rs})
#     print(r""">>>
# The homogenized elastic coefficient $E_{ijkl}$, symmetric storage
# with rows, columns in 11, 22, 12 ordering:""")
    #print(c_e)
    return sff.youngpoisson_from_stiffness(c_e,plane='strain')

if __name__ == '__main__':
    meshes = ["Input/cercle", "Input/ellipse", "Input/rectangle"]
    lc = [1.5e-2]
    res = {}
    l = 2
    for mesh in meshes:
        key = mesh[6:]
        res[key] = []
        for i in range(l):
            os.rename(mesh + str(lc[0]) + 'n' + str(i) + '.h5', 'Input/meshfile.h5')
            res[key].append(main())
            os.rename('Input/meshfile.h5', mesh + str(lc[0]) + 'n' + str(i) + '.h5')
    print(res)