#!/usr/bin/env python
# -*- Python -*-
#
# Copyright 2003 Free Software Foundation, Inc.
# 
# This file is part of GNU Radio
# 
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
# 
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 


# Edit the gpif.c file generated by the Cypress GPIF Designer Tool and
# produce usrp_gpif.c, and usrp_gpif_inline.h, files suitable for our
# uses.

import re
import string
import sys

def check_flow_state (line, flow_state_dict):
    mo = re.match (r'/\* Wave (\d) FlowStates \*/ (.*),', line)
    if mo:
        wave = int (mo.group (1))
        data = mo.group (2)
        split = data.split (',', 8)
        v = map (lambda x : int (x, 16), split)
        # print "%s, %s" % (wave, data)
        # print "split: ", split
        # print "v    : ", v
        flow_state_dict[wave] = v


def delta (xseq, yseq):
    # set subtraction
    z = []
    for x in xseq:
        if x not in yseq:
            z.append (x)
    return z
    

def write_define (output, name, pairs):
    output.write ('#define %s()\t\\\n' % name)
    output.write ('do {\t\t\t\t\t\\\n')
    for reg, val in pairs:
        output.write ('%14s = 0x%02x;\t\t\t\\\n' % (reg, val))
    output.write ('} while (0)\n\n')
    
def write_inlines (output, dict):
    regs = ['FLOWSTATE', 'FLOWLOGIC', 'FLOWEQ0CTL', 'FLOWEQ1CTL', 'FLOWHOLDOFF',
            'FLOWSTB', 'FLOWSTBEDGE', 'FLOWSTBHPERIOD', 'GPIFHOLDAMOUNT']

    READ_CTRL_FLOW_STATE = 0
    WRITE_CTRL_FLOW_STATE = 1
    READ_DATA_FLOW_STATE = 2
    WRITE_DATA_FLOW_STATE = 3

    read_data_info = zip (regs, dict[READ_DATA_FLOW_STATE])
    write_data_info = zip (regs, dict[WRITE_DATA_FLOW_STATE])
    read_ctrl_info = zip (regs, dict[READ_CTRL_FLOW_STATE])
    write_ctrl_info = zip (regs, dict[WRITE_CTRL_FLOW_STATE])
    
    output.write ('''/*
 * Machine generated by "edit-gpif".  Do not edit by hand.
 */

''')
    write_define (output, 'setup_flowstate_common', read_data_info) #assumes that the same registers will change, this isn't really good
    write_define (output, 'setup_flowstate_data_read', delta (read_data_info, write_data_info))
    write_define (output, 'setup_flowstate_data_write', delta (write_data_info, read_data_info))
    write_define (output, 'setup_flowstate_ctrl_read', delta (read_ctrl_info, write_ctrl_info))
    write_define (output, 'setup_flowstate_ctrl_write', delta (write_ctrl_info, read_ctrl_info))
    

def edit_gpif (input_name, output_name, inline_name):
    input = open (input_name, 'r')
    output = open (output_name, 'w')
    inline = open (inline_name, 'w')
    flow_state_dict = {}

    output.write ('''/*
 * Machine generated by "edit-gpif".  Do not edit by hand.
 */

''')
    
    while 1:
        line = input.readline ()
        line = string.replace (line, '\r','')
        line = re.sub (r' *$', r'', line)

        check_flow_state (line, flow_state_dict)

        line = re.sub (r'#include', r'// #include', line)
        line = re.sub (r'xdata ', r'', line)
        if re.search (r'GpifInit', line):
            break
        
        output.write (line)

    output.close ()
    write_inlines (inline, flow_state_dict)
    inline.close ()


# gpif.c usrp_gpif.c usrp_gpif_inline.h
edit_gpif (sys.argv[1], sys.argv[2], sys.argv[3])
