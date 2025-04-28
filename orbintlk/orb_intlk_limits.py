# This is a python3 file:
#     - The object defined here is a python3 dictionary.
#     - It can be imported into python3 code with the code:
#         with open('/pathto/orb_intlk_limits.py') as txt:
#             return_dict = ast.literal_eval(txt.read())
#     where `ast` is a python3 module available in `pip3`
#
# This file summarizes the orbit interlock limits for position and angle for
# all relevant IDs currently in operation at SIRIUS.

{  # units in [um] and [urad]  # noqa: B018
    '06SB': {'psx': 600, 'psy': 600, 'agx': 300, 'agy': 300},  # VPU
    '08SB': {'psx': 1000, 'psy': 1000, 'agx': 400, 'agy': 200},  # IVU's
    '10SB': {'psx': 500, 'psy': 500, 'agx': 400, 'agy': 400},  # Delta
    '14SB': {'psx': 1000, 'psy': 1000, 'agx': 400, 'agy': 200},  # IVU's
}
