from pymol import cmd

def setdefaults():
    cmd.color('black', 'name C')
    cmd.color('yellow', 'name S')
    cmd.hide('lines', 'all')
    cmd.set('stick_radius',0.1)
    cmd.set('sphere_scale',0.25)
    cmd.set('sphere_scale',1,'name au')
    cmd.set('sphere_scale',1,'name ag')
    cmd.color('grey90', 'name ag')
    cmd.color('gold', 'name au')
    cmd.set('antialias',2)
    cmd.set('opaque_background',0)
    cmd.set('ray_trace_gain',5)
    cmd.set('depth_cue',0)
    cmd.show('sticks','all')
    cmd.show('spheres','all')
cmd.extend('setdefaults',setdefaults)
