# This is a python3 file:
#     - The object defined here is a python3 dictionary.
#     - It can be imported into python3 code with the code:
#         `return_dict = ast.literal_eval('/pathto/high-level-triggers.py')`
#     where `ast` is a python3 module available in `pip3`
#
# This file define the triggers used in SIRIUS timing test banch setup.
#
#
# To run the high level timing windows on your PC you must have:
#  - sirius mamba environment installed
#  - a SIRIUS Control System constants runing in your local machine or some
#    other machine you have access.
#
# If the requirements above are satisfied, you should then change the branch
# of the active control-system-constants repo to the branch:
#    git checktout timesys-testbench-configs
#
# 1-) Then, on a clean terminal, you must activate the sirius enviroment
#    mamba activate sirius
#
# 2-) Export the enviroment variables defining where to find the low level PVs:
#    export EPICS_CA_ADDR_LIST="$EPICS_CA_ADDR_LIST 10.10.10.126 10.20.26.78 10.20.26.82"
#    export EPICS_PVA_ADDR_LIST="$EPICS_PVA_ADDR_LIST 10.10.10.126 10.20.26.78 10.20.26.82"
#
# 3-) Export the environment variable pointing to the right cs-consts server:
#    export SIRIUS_URL_CONSTS=http://localhost
#
# 4-) Run the high level testing IOC on you machine:
#    sirius-ioc-as-ti-control.py
#
# 5-) Open a new terminal and run steps 1-3 again. Then run the windows:
#    sirius-hla-as-ti-control.py
{
    'AS-01:TI-TestEVR1': {
        'database': {
            'Src': {'enums': ('InjSI', )},
            'Delay': {'value': 0.0},
            'RFDelayType': {'value': 1, 'states': (2, 0)},
            'NrPulses': {'value': 1, 'hilim': 2, 'high': 2, 'hihi': 2},
            'Duration': {'value': 150, 'hilim': 550, 'high': 550, 'hihi': 600},
            'State': {'value': 0},
            'Polarity': {'value': 0, 'states': (0, 2)},
            },
        'channels': (
            'DE-23RaBPM:DU-DummyEVR1:IN',
            'DE-23RaBPM:DU-DummyEVR2:IN',
            ),
        },
    'AS-01:TI-TestEVR2': {
        'database': {
            'Src': {'enums': ('InjSI', )},
            'Delay': {'value': 0.0},
            'RFDelayType': {'value': 1, 'states': (2, 0)},
            'NrPulses': {'value': 1, 'hilim': 2, 'high': 2, 'hihi': 2},
            'Duration': {'value': 150, 'hilim': 550, 'high': 550, 'hihi': 600},
            'State': {'value': 0},
            'Polarity': {'value': 1, 'states': (2, 0)},
            },
        'channels': (
            'DE-23RaBPM:DU-DummyEVR3:IN',
            'DE-23RaBPM:DU-DummyEVR4:IN',
            ),
        },
    'AS-01:TI-TestEVE1': {
        'database': {
            'Src': {'enums': ('InjSI', )},
            'Delay': {'value': 0.0},
            'RFDelayType': {'value': 1, 'states': (2, 0)},
            'NrPulses': {'value': 1, 'hilim': 2, 'high': 2, 'hihi': 2},
            'Duration': {'value': 150, 'hilim': 550, 'high': 550, 'hihi': 600},
            'State': {'value': 0},
            'Polarity': {'value': 0, 'states': (0, 2)},
            },
        'channels': (
            'DE-23RaBPM:DU-DummyEVE1:IN',
            'DE-23RaBPM:DU-DummyEVE2:IN',
            ),
        },
    'AS-01:TI-TestEVE2': {
        'database': {
            'Src': {'enums': ('InjSI', )},
            'Delay': {'value': 0.0},
            'RFDelayType': {'value': 1, 'states': (2, 0)},
            'NrPulses': {'value': 1, 'hilim': 2, 'high': 2, 'hihi': 2},
            'Duration': {'value': 150, 'hilim': 550, 'high': 550, 'hihi': 600},
            'State': {'value': 0},
            'Polarity': {'value': 1, 'states': (2, 0)},
            },
        'channels': (
            'DE-23RaBPM:DU-DummyEVE3:IN',
            'DE-23RaBPM:DU-DummyEVE4:IN',
            ),
        },
    'AS-01:TI-TestAMC1': {
        'database': {
            'Src': {'enums': ('InjSI', )},
            'Delay': {'value': 0.0},
            'RFDelayType': {'value': 1, 'states': (2, 0)},
            'NrPulses': {'value': 1, 'hilim': 2, 'high': 2, 'hihi': 2},
            'Duration': {'value': 150, 'hilim': 550, 'high': 550, 'hihi': 600},
            'State': {'value': 0},
            'Polarity': {'value': 1, 'states': (2, 0)},
            },
        'channels': (
            'DE-23RaBPM:DU-DummyAMC3:IN',
            'DE-23RaBPM:DU-DummyAMC4:IN',
            'DE-23RaBPM:DU-DummyAMC5:IN',
            'DE-23RaBPM:DU-DummyAMC6:IN',
            ),
        },
    'AS-01:TI-TestAMC2': {
        'database': {
            'Src': {'enums': ('InjSI', )},
            'Delay': {'value': 0.0},
            'RFDelayType': {'value': 1, 'states': (2, 0)},
            'NrPulses': {'value': 1, 'hilim': 2, 'high': 2, 'hihi': 2},
            'Duration': {'value': 150, 'hilim': 550, 'high': 550, 'hihi': 600},
            'State': {'value': 0},
            'Polarity': {'value': 1, 'states': (2, 0)},
            },
        'channels': (
            'DE-23RaBPM:DU-DummyAMC1:IN',
            'DE-23RaBPM:DU-DummyAMC2:IN',
            ),
        },
}
