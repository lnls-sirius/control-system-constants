# This is a python3 file:
#     - The object defined here is a python3 dictionary.
#     - It can be imported into python3 code with the code:
#         with open('/pathto/high-level-events.py') as txt:
#             return_dict = ast.literal_eval(txt.read())
#     where `ast` is a python3 module available in `pip3`
#
# This file maps the aliases of the events used in Sirius
# operation system with their code name.
{  # noqa: B018
    'Linac': 'Evt01', 'InjBO': 'Evt02',
    'InjSI': 'Evt03', 'RmpBO': 'Evt04',
    'Cycle': 'Evt05', 'Study': 'Evt06',
    'OrbSI': 'Evt07', 'CplSI': 'Evt08',
    'TunSI': 'Evt09', 'FOFBS': 'Evt10',

    'Dsbld':  'Evt00',
    'Intlk': 'Evt117', 'ItlkR': 'Evt118',
    'LLRFA': 'Evt121', 'LLRFB': 'Evt122',
    'RFKll': 'Evt124',
    'DCT13': 'Evt125', 'DCT14': 'Evt126',
    'PsMtm': 'Evt132',
}
