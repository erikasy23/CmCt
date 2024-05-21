#!/usr/bin/env python

import ipywidgets as widgets
from ipywidgets import VBox, HBox, Text, HTML, Dropdown, Button

def CMCT_gaivimetry_select_data_parameters():

    # User Input Widgets
    loginname_widget = Text(value='rbasnet', description='Login Name:')

    # Collect all widgets in a vertical layout
    all_widgets = VBox([loginname_widget])

    # Display all widgets
    display(all_widgets)
    
    # Return a dictionary of the current values
    return {     
        'loginname_widget': loginname_widget.value
    }







