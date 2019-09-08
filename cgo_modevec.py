'''
modified http://pymolwiki.org/index.php/cgo_arrow
to work for modevectors from the chem package

(c) 2013 Thomas Holder, Schrodinger Inc.

License: BSD-2-Clause
modified by: Jeff Becca
'''

from pymol import cmd, cgo, CmdException


def cgo_modevec(atom1='pk1', atom2='pk2', radius=0.05, gap=0.0, hlength=-1, hradius=-1,
              color='green', name='',scalefactor=10.0, cutoff=0.6, transparency=1.0): #was 0.6cut 12 scale
## scalefactor was 10.0 
    '''
DESCRIPTION

    Create a CGO mode vector starting at atom1 and pointint in atom2 displacement

ARGUMENTS

    atom1 = string: single atom selection or list of 3 floats {default: pk1}

    atom2 = string: displacement to atom1 for modevec

    radius = float: arrow radius {default: 0.5}

    gap = float: gap between arrow tips and the two atoms {default: 0.0}

    hlength = float: length of head

    hradius = float: radius of head

    color = string: one or two color names {default: blue red}

    name = string: name of CGO object
    
    scalefactor = scale how big of an arrow to make. Default 5

    transparency = 0.0 ~ 1.0, default=1.0 means being totally opaque
    '''
    from chempy import cpv

    radius, gap = float(radius), float(gap)
    hlength, hradius = float(hlength), float(hradius)
    scalefactor, cutoff = float(scalefactor), float(cutoff)
    transparency = float(transparency)
    try:
        color1, color2 = color.split()
    except:
        color1 = color2 = color
    color1 = list(cmd.get_color_tuple(color1))
    color2 = list(cmd.get_color_tuple(color2))

    def get_coord(v):
        if not isinstance(v, str):
            return v
        if v.startswith('['):
            return cmd.safe_list_eval(v)
        return cmd.get_atom_coords(v)

    xyz1 = get_coord(atom1)
    xyz2 = get_coord(atom2)
    newxyz2 = cpv.scale(xyz2, scalefactor)
    newxyz2 = cpv.add(newxyz2, xyz1)
    xyz2 = newxyz2
#    xyz2 = xyz2[0]*scalefactor, xyz2[1]*scalefactor, xyz2[2]*scalefactor
    normal = cpv.normalize(cpv.sub(xyz1, xyz2))

    if hlength < 0:
        hlength = radius * 3.0
    if hradius < 0:
        hradius = hlength * 0.6

    if gap:
        diff = cpv.scale(normal, gap)
        xyz1 = cpv.sub(xyz1, diff)
        xyz2 = cpv.add(xyz2, diff)

    xyz3 = cpv.add(cpv.scale(normal, hlength), xyz2)

# dont draw arrow if distance is too small
    distance = cpv.distance(xyz1, xyz2)
    if distance <= cutoff:
        return

#### generate transparent arrows; 
#### The original codes are the next block
####  --Ran
    obj = [25.0, transparency, 9.0] + xyz1 + xyz3 + [radius] + color1 + color2 + \
          [25.0, transparency, 27.0] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
          [1.0, 0.0]

#    obj = [cgo.CYLINDER] + xyz1 + xyz3 + [radius] + color1 + color2 + \
#          [cgo.CONE] + xyz3 + xyz2 + [hradius, 0.0] + color2 + color2 + \
#          [1.0, 0.0]

    if not name:
        name = cmd.get_unused_name('arrow')

    cmd.load_cgo(obj, name)

cmd.extend('cgo_modevec', cgo_modevec)
