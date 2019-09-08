from pymol.cgo import *
from pymol import cmd
from pymol.vfont import plain

def betteraxes():
    w = 0.06 # cylinder width 
    l = 0.75 # cylinder length
    h = 0.25 # cone hight
    d = w * 1.618 # cone base diameter
    
    obj = [CYLINDER, 0.0, 0.0, 0.0,   l, 0.0, 0.0, w, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0,   l, 0.0, w, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0,
           CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   l, w, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
           CONE,   l, 0.0, 0.0, h+l, 0.0, 0.0, d, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 
           CONE, 0.0,   l, 0.0, 0.0, h+l, 0.0, d, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 
           CONE, 0.0, 0.0,   l, 0.0, 0.0, h+l, d, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0]

    cyl_text(obj,plain,[1.1,-0.04,0.],'X',0.025,axes=[[0.125,0,0],[0,0.125,0],[0,0,0.125]],color=[1.0,0.0,0.0])
    cyl_text(obj,plain,[-0.04,1.1,0.],'Y',0.025,axes=[[0.125,0,0],[0,0.125,0],[0,0,0.125]], color=[0.0, 1.0, 0.0])
    cyl_text(obj,plain,[0.,0.,1.1],'Z',0.025,axes=[[0.125,0,0],[0,0.125,0],[0,0,0.125]], color=[0.0, 0.0, 1.0])

    cmd.load_cgo(obj, 'axes')

cmd.extend('betteraxes',betteraxes)
