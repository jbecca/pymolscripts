'''
Calculation of intertia tensor taken from
http://www.pymolwiki.org/index.php/inertia_tensor
 
(c) August 2010 by Mateusz Maciejewski
matt (at) mattmaciejewski . com

License: MIT


Rest of code translates coordinates to yeild a COM of [0,0,0] 
and then aligns one of the principle axes of rotation to the Z unit vector.
'''

from pymol.cgo import *
from pymol import cmd


def alignZ(selection, name="tensor", state=1, scaling=0, quiet=1):
    """
DESCRIPTION

    This script will draw the inertia tensor of the selection.

ARGUMENTS

    selection = string: selection for the atoms included in the tensor calculation

    name = string: name of the tensor object to be created {default: "tensor"}

    state = int: state/model in the molecule object used in the tensor calculation

    scaling = int {0, 1, or 2}: 0 for no scaling of the inertia axes, 1 for scaling
    according to the molecular shape, 2 for scaling according to the eigenvalues 
    {default: 0}

EXAMPLE

    PyMOL> run inertia_tensor.py
    PyMOL> tensor molecule_object & i. 2-58+63-120 & n. n+ca+c, "tensor_model5_dom2", 5, 1

NOTES

    Requires numpy.
    """

    import numpy

    totmass = 0.0
    x_com, y_com, z_com = 0, 0, 0

    model = cmd.get_model(selection, state)

    for a in model.atom:

        x_com += a.coord[0] * a.get_mass()
        y_com += a.coord[1] * a.get_mass()
        z_com += a.coord[2] * a.get_mass()
        totmass += a.get_mass()

    x_com /= totmass
    y_com /= totmass
    z_com /= totmass

    if not int(quiet):
        print()
        print("Center of mass: ")
        print()
        print(x_com, y_com, z_com)
    cmd.translate([-x_com, -y_com, -z_com], selection)

    # Go back and rewrite the COM section to be a function, then it can
    # be called again to ensure the translation was correct
    x_com = 0.0
    y_com = 0.0
    z_com = 0.0

    I = []

    for index in range(9):
        I.append(0)

    for a in model.atom:

        temp_x, temp_y, temp_z = a.coord[0], a.coord[1], a.coord[2]
        temp_x -= x_com
        temp_y -= y_com
        temp_z -= z_com

        I[0] += a.get_mass() * (temp_y ** 2 + temp_z ** 2)
        I[4] += a.get_mass() * (temp_x ** 2 + temp_z ** 2)
        I[8] += a.get_mass() * (temp_x ** 2 + temp_y ** 2)
        I[1] -= a.get_mass() * temp_x * temp_y
        I[3] -= a.get_mass() * temp_x * temp_y
        I[2] -= a.get_mass() * temp_x * temp_z
        I[6] -= a.get_mass() * temp_x * temp_z
        I[5] -= a.get_mass() * temp_y * temp_z
        I[7] -= a.get_mass() * temp_y * temp_z

    tensor = numpy.array([(I[0:3]), (I[3:6]), (I[6:9])])
    vals, vects = numpy.linalg.eig(tensor)  # they come out unsorted, so the command below is needed

    eig_ord = numpy.argsort(vals)  # a thing to note is that here COLUMN i corrensponds to eigenvalue i.

    ord_vals = vals[eig_ord]
    ord_vects = vects[:, eig_ord].T

    if not int(quiet):
        print()
        print("Inertia tensor z, y, x eigenvalues:")
        print()
        print(ord_vals)
        print()
        print("Inertia tensor z, y, x eigenvectors:")
        print()
        print(ord_vects)

    if int(scaling) == 0:
        norm_vals = [sum(numpy.sqrt(ord_vals / totmass)) / 3 for i in range(3)]

    elif int(scaling) == 1:
        normalizer = numpy.sqrt(max(ord_vals) / totmass)
        norm_vals = normalizer / numpy.sqrt(ord_vals / totmass) * normalizer
        norm_vals = norm_vals / (max(norm_vals) / min(norm_vals))

    elif int(scaling) == 2:
        normalizer = numpy.sqrt(max(ord_vals) / totmass)
        norm_vals = numpy.sqrt(ord_vals / totmass)

    start = [x_com, y_com, z_com]
    ends = [[(norm_vals[0] - 1) * ord_vects[0][0], (norm_vals[0] - 1) * ord_vects[0][1], (norm_vals[0] - 1) * ord_vects[0][2]],
          [(norm_vals[1] - 1) * ord_vects[1][0], (norm_vals[1] - 1) * ord_vects[1][1], (norm_vals[1] - 1) * ord_vects[1][2]],
          [(norm_vals[2] - 1) * ord_vects[2][0], (norm_vals[2] - 1) * ord_vects[2][1], (norm_vals[2] - 1) * ord_vects[2][2]]]

    import numpy as np 
    zvec = np.array([[0],[0],[1]])
    newvec = np.array([[(norm_vals[0]-1) * ord_vects[0][2]],
                       [(norm_vals[1]-1) * ord_vects[1][2]],
                       [(norm_vals[2]-1) * ord_vects[2][2]]])
    print(newvec)
    norm = np.linalg.norm(newvec)
    newvecnorm = np.divide(newvec, norm)
    newvec = newvecnorm
    print(newvec)

    #create the rotation matrix, cross product, and dot prduct
    rotmat = np.zeros((3,3))
    cross = np.cross(newvec, zvec, axis=0)
    crossnorm = np.linalg.norm(cross)
    dot = np.dot(np.transpose(newvec), zvec)

    vxmat = np.zeros((3,3))
    vxmat[0,1] = -cross[2]
    vxmat[0,2] = cross[1]
    vxmat[1,0] = cross[2]
    vxmat[1,2] = -cross[0]
    vxmat[2,0] = -cross[1]
    vxmat[2,1] = cross[0]

    thirdterm = np.matmul(vxmat, vxmat) * (1/(1+dot))

    rotmat = np.diag((1,1,1)) + vxmat + thirdterm
    newvectrans = np.matmul(rotmat, newvec)
    print("New vector after alignment: %s" % newvectrans)
    cmd.transform_selection(selection, ([rotmat[0,0], rotmat[0,1], rotmat[0,2], 0.0, 
                                        rotmat[1,0], rotmat[1,1], rotmat[1,2], 0.0,
                                        rotmat[2,0], rotmat[2,1], rotmat[2,2], 0.0,
                                        0.0,            0.0,        0.0,        1.0]))

cmd.extend("aligntoZ", alignZ)