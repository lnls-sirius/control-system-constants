# This is a python3 file:
#     - The object defined here is a python3 dictionary.
#     - It can be imported into python3 code with the code:
#         `return_dict = ast.literal_eval('/pathto/high-level-events.py')`
#     where `ast` is a python3 module available in `pip3`
#
# This file summarizes the orbit interlock limits for position and angle for
# all relevant IDs currently in operation at SIRIUS.

{  # units in [um] and [urad]
    'default': dict(psx=1500, psy=1500, agx=500, agy=500),  # Rest of the ring
    '06SB': dict(psx=600, psy=600, agx=300, agy=300),  # VPU
    '08SB': dict(psx=1000, psy=1000, agx=400, agy=200),  # IVU's
    '10SB': dict(psx=500, psy=500, agx=400, agy=400),  # Delta
    '14SB': dict(psx=1000, psy=1000, agx=400, agy=200),  # IVU's
}
