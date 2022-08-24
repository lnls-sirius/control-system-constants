# This is a python3 file:
#     - The object defined here is a python3 dictionary.
#     - It can be imported into python3 code with the code:
#         `return_dict = ast.literal_eval('/pathto/crate-resource-meaning.py')`
#     where `ast` is a python3 module available in `pip3`
#
# This file maps the meaning of the crate backplane resources.
{
    'MTCA.4 backplane resource': [
        'M-LVDS Trigger 1', 'M-LVDS Trigger 2', 'M-LVDS Trigger 3',
        'M-LVDS Trigger 4', 'M-LVDS Trigger 5', 'M-LVDS Trigger 6',
        'M-LVDS Trigger 7', 'M-LVDS Trigger 8',
        'Clk1/MCH1 (TCLKA)', 'Clk2/MCH1 (TCLKB)',
    ],
    'Channel': [
        'CH0', 'CH1', 'CH2', 'CH3', 'CH4', 'CH5', 'CH6', 'CH7', '', '',
    ],
    'Type': [
        'Trigger', 'Trigger', 'Trigger', 'Trigger', 'Trigger', 'Trigger',
        '(reserved)', 'Clock', 'Clock', '(reserved)',
    ],
    'Description': [
        'Machine 1', 'Machine 2', 'Beam Loss (Timing Receiver to BPM)',
        'Beam Loss (BPM to Timing Receiver)', 'Orbit Interlock',
        'FOFB DCC Sync', '(reserved)', 'SOFB (Monit) clock',
        'BPM reference clock', '(reserved)',
    ],
    'Source': [
        'Timing Receiver', 'Timing Receiver', 'Timing Receiver', 'All BPMs',
        'All BPMs', 'Timing Receiver', 'Timing Receiver', 'Timing Receiver',
        'Timing Receiver (TCLKA)', '(reserved)',
    ],
    'Details': [
        '', '', '', '', '', '', '', 'Freq. RF/19872000', 'Freq. 5/36*RF',
        '(reserved)',
    ]
}
