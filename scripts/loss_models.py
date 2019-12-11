import pypsa
import numpy as np
import subprocess
from tempfile import mkstemp
#from vresutils import load as vload

import pandas as pd
import scipy.io
import re
import logging
from pypsa.opt import LConstraint, l_constraint, LExpression
from six import iteritems, itervalues, string_types
import pyomo.environ as penv
from pyomo.environ import Constraint, Objective, Var, ComponentUID
import pandas as pd


def post_cosine(network, snapshots,duals):
    nrOfIntervalls = 10
    passive_branches = network.passive_branches()
    loss = pd.DataFrame(0, index=snapshots, columns = list((i[1] for i in network.passive_branches().index))) 
    for branch in passive_branches.index:
        sub = passive_branches.at[branch,"sub_network"]
        attribute_r = "r_pu_eff" if network.sub_networks.at[sub,"carrier"] == "DC" else "x_pu_eff"
        attribute_s = "s_nom"
        conductance = 1/passive_branches.at[branch, attribute_r]
        bt = branch[0]
        bn = branch[1]
        xU = (passive_branches.at[branch,attribute_s]+6800)/passive_branches.at[branch, attribute_r] + 6800
        for sn in snapshots:
            for i in list(range(nrOfIntervalls)): 
                lower = ((i/nrOfIntervalls) * xU)
                max_val = lower+(xU/nrOfIntervalls)
                slope = 2*conductance*np.sin(lower)
                loss.loc[sn,bt]=((1/2)*slope*(network.model.delta_angle_positive[bt,bn,sn,i]+network.model.delta_angle_negative[bt,bn,sn,i]))
                loss.loc[sn,bn]=((1/2)*slope*(network.model.delta_angle_positive[bt,bn,sn,i]+network.model.delta_angle_negative[bt,bn,sn,i]))    
    network.lines_t["loss"] =  loss
def cosine(network,snapshots):
    passive_branches = network.passive_branches()
    nrOfIntervalls = 10
    network.model.delta_angle_positive = Var(list(passive_branches.index), snapshots,list(range(nrOfIntervalls)))
    network.model.delta_angle_negative = Var(list(passive_branches.index), snapshots,list(range(nrOfIntervalls)))
    network.model.difference_constraint = penv.ConstraintList()
    intervalls_upper = {}
    intervalls_lower = {}
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch, "bus0"]
        bus1 = passive_branches.at[branch, "bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]
        attribute_r = "r_pu_eff" if network.sub_networks.at[sub,"carrier"] == "AC" else "x_pu_eff"
        attribute_s = "s_nom"
        conductance = 1/passive_branches.at[branch, attribute_r]
        xU = (passive_branches.at[branch,attribute_s]+6800)/passive_branches.at[branch, attribute_r] + 6800
        bt = branch[0]
        bn = branch[1]
        for sn in snapshots:
            for i in list(range(nrOfIntervalls)): 
                lower = ((i/nrOfIntervalls) * xU)
                max_val = lower+(xU/nrOfIntervalls)
                slope = 2*conductance*np.sin(lower)
                intervalls_lower[bt,bn,sn,i] = (LConstraint(LExpression([(1,network.model.delta_angle_positive[bt,bn,sn,i])]), ">=", LExpression()))
                intervalls_upper[bt,bn,sn,i] = (LConstraint(LExpression([(1,network.model.delta_angle_negative[bt,bn,sn,i])],(-max_val)), "<=", LExpression()))  
            losses = sum(slope*(network.model.delta_angle_positive[bt,bn,sn,i]+network.model.delta_angle_negative[bt,bn,sn,i]) for i in list(range(nrOfIntervalls)))
            network.model.power_balance[bus0,sn]._body -=((1/2)*losses)
            network.model.power_balance[bus1,sn]._body -=((1/2)*losses)    
            network.model.difference_constraint.add((network.model.voltage_angles[bus0,sn] - network.model.voltage_angles[bus1,sn]) ==  (sum((network.model.delta_angle_positive[bt,bn,sn,i] - network.model.delta_angle_negative[bt,bn,sn,i]) for i in list(range(nrOfIntervalls))))      )
    l_constraint(network.model,"intervalls_lower",intervalls_lower,list(passive_branches.index),snapshots, list(range(nrOfIntervalls)))
    l_constraint(network.model,"intervalls_upper",intervalls_upper,list(passive_branches.index), snapshots, list(range(nrOfIntervalls)))   
def post_quadratic(network, snapshots,duals):
    
    passive_branches = network.passive_branches()
    loss = pd.DataFrame(0, index=snapshot, columns=list((i[1] for i in network.passive_branches().index))) 
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch, "bus0"]
        bus1 = passive_branches.at[branch, "bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        attribute_r = "r_pu_eff" if network.sub_networks.at[sub,"carrier"] == "DC" else "x_pu_eff"
        conductance = 1/passive_branches.at[branch, attribute_r]
        bt = branch[0]
        bn = branch[1]
        for sn in snapshots: 
            loss.loc[sn,bn] = (conductance*network.model.delta_angle_sq[bt,bn,sn])
    network.lines_t["loss"] =  loss
def quadratic(network,snapshots):
    nrOfIntervalls = 10
    passive_branches = network.passive_branches()
    network.model.delta_angle = Var(list(passive_branches.index), snapshots)
    network.model.delta_angle_sq = Var(list(passive_branches.index), snapshots)  
    differences = {}
    envelope = {}
    envelope_index = ["wov1","xup","xlow"]
    intervalls = {}
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch, "bus0"]
        bus1 = passive_branches.at[branch, "bus1"]
        bt = branch[0]
        bn = branch[1]
        sub = passive_branches.at[branch,"sub_network"]
        attribute_s = "s_nom"
        attribute_r = "r_pu_eff" if network.sub_networks.at[sub,"carrier"] == "AC" else "x_pu_eff"
        conductance = 1/passive_branches.at[branch, attribute_r]
        r = passive_branches.at[branch, attribute_r]
        xU = (passive_branches.at[branch,attribute_s]+6800)/passive_branches.at[branch, attribute_r]
        xL = -xU
        for sn in snapshots:
            lhs = LExpression([(1,network.model.voltage_angles[bus0,sn]),(-1,network.model.voltage_angles[bus1,sn]),(-1,network.model.delta_angle[bt,bn,sn])])
            differences[bt,bn,sn] = LConstraint(lhs,"==",LExpression())
            #bounds 
            envelope[bt,bn,sn,"xup"]=(LConstraint(LExpression([(1,network.model.delta_angle[bt,bn,sn])],(-xU)), "<=", LExpression()))
            envelope[bt,bn,sn, "xlow"]=(LConstraint(LExpression([(1,network.model.delta_angle[bt,bn,sn])],(xU)),">=", LExpression()))                                                       
            #approximation from above
            envelope[bt,bn,sn,"wov1"]=(LConstraint(LExpression([(1,network.model.delta_angle_sq[bt,bn,sn])]), "<=", LExpression([(xU,network.model.delta_angle[bt,bn,sn]),(xL,network.model.delta_angle[bt,bn,sn])],(-xU*xL))))
            #appoximation from below in  intervalls
            for i in list(range(nrOfIntervalls+1)):
                lower = ((i/nrOfIntervalls) * xU)
                intervalls[i,bt,bn,sn] = (LConstraint(LExpression([(1,network.model.delta_angle_sq[bt,bn,sn])]),">=",LExpression([(2*lower,network.model.delta_angle[bt,bn,sn])],(-lower*lower))))
            #add p_sq to nodal power balance
            #use of ._body because of pyomo bug
            network.model.power_balance[bus0,sn]._body -=((1/2)*conductance*network.model.delta_angle_sq[bt,bn,sn])
            network.model.power_balance[bus1,sn]._body -=((1/2)*conductance*network.model.delta_angle_sq[bt,bn,sn])
    l_constraint(network.model,"mccormick_envelope_delta_angle", envelope, list(passive_branches.index),snapshots,envelope_index)
    l_constraint(network.model,"differences", differences, list(passive_branches.index),snapshots)
    l_constraint(network.model,"envelope_intervalls",intervalls,list(range(nrOfIntervalls+1)), list(passive_branches.index), snapshots)    
def post_lldc(network,snapshots,duals):
    passive_branches = network.passive_branches()
    loss = pd.DataFrame(0, index=snapshots, columns=list((i[1] for i in passive_branches.index)))
    for branch in passive_branches.index:
        sub = passive_branches.at[branch,"sub_network"]
        attribute_r = "r_pu_eff" if network.sub_networks.at[sub,"carrier"] == "DC" else "x_pu_eff"
        r = passive_branches.at[branch, attribute_r]
        bt = branch[0]
        bn = branch[1]
        for sn in snapshots: 
            loss.loc[sn,bn] = (r*network.model.passive_branch_p_sq_out[bt,bn,sn].value) + (r*network.model.passive_branch_p_sq_in[bt,bn,sn].value)
    network.lines_t["loss"] =  loss
        
#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6345342
def lldc(network,snapshots):
    nrOfintervalls = 10
    passive_branches = network.passive_branches()
    network.model.passive_branch_p_out = Var(list(passive_branches.index), snapshots, bounds = (0,float('inf')))
    network.model.passive_branch_p_in = Var(list(passive_branches.index), snapshots, bounds = (0,float('inf')))
    network.model.passive_branch_p_sq_in = Var(list(passive_branches.index), snapshots, bounds = (0,float('inf')))
    network.model.passive_branch_p_sq_out = Var(list(passive_branches.index), snapshots, bounds = (0,float('inf')))
    
    losses = {}
    approx_index = ["in","out"]
    envelope = {}
    envelope_index = ["xup","wov1"]
    intervalls = {}
    
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch, "bus0"]
        bus1 = passive_branches.at[branch, "bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        attribute_r = "r_pu_eff" if network.sub_networks.at[sub,"carrier"] == "AC" else "x_pu_eff"
        if(passive_branches.at[branch,"s_nom_extendable"]):
            xU = passive_branches.at[branch,"s_nom"] +6800
        else:
            xU = passive_branches.at[branch,"s_nom"]
        r = passive_branches.at[branch, attribute_r]
        bt = branch[0]
        bn = branch[1]
        xL = 0
        for sn in snapshots:
            lhs = LExpression([(1,network.model.passive_branch_p_out[bt,bn,sn]),(-1,network.model.passive_branch_p_in[bt,bn,sn]),(-1,network.model.passive_branch_p[bt,bn,sn])])
            losses[bt,bn,sn] = LConstraint(lhs,"==",LExpression())
            #bounds
            envelope["in",bt,bn,sn,"xup"]=(LConstraint(LExpression([(1,network.model.passive_branch_p_in[bt,bn,sn]),(r,network.model.passive_branch_p_sq_in[bt,bn,sn]),(r,network.model.passive_branch_p_sq_out[bt,bn,sn])],(-xU)), "<=", LExpression()))
            envelope["out",bt,bn,sn,"xup"]=(LConstraint(LExpression([(1,network.model.passive_branch_p_out[bt,bn,sn]),(r,network.model.passive_branch_p_sq_in[bt,bn,sn]),(r,network.model.passive_branch_p_sq_out[bt,bn,sn])],(-xU)), "<=", LExpression()))
            #approximation from above      
            envelope["in",bt,bn,sn,"wov1"]=(LConstraint(LExpression([(1,network.model.passive_branch_p_sq_in[bt,bn,sn])]), "<=", LExpression([(xU,network.model.passive_branch_p_in[bt,bn,sn]),(xL,network.model.passive_branch_p_in[bt,bn,sn])],(-xU*xL))))
            envelope["out",bt,bn,sn,"wov1"]=(LConstraint(LExpression([(1,network.model.passive_branch_p_sq_out[bt,bn,sn])]), "<=", LExpression([(xU,network.model.passive_branch_p_out[bt,bn,sn]),(xL,network.model.passive_branch_p_out[bt,bn,sn])],(-xU*xL))))
            #appoximation from below in 10 intervalls
            for i in list(range(nrOfintervalls+1)):
                lower = ((i/nrOfintervalls) * xU)
                intervalls["in",i,bt,bn,sn] = (LConstraint(LExpression([(1,network.model.passive_branch_p_sq_in[bt,bn,sn])]),">=",LExpression([(2*lower,network.model.passive_branch_p_in[bt,bn,sn])],(-lower*lower))))
                intervalls["out",i,bt,bn,sn] = (LConstraint(LExpression([(1,network.model.passive_branch_p_sq_out[bt,bn,sn])]),">=",LExpression([(2*lower,network.model.passive_branch_p_out[bt,bn,sn])],(-lower*lower))))
            #add p_sq to nodal power balance
            #use of ._body because of pyomo bug
            network.model.power_balance[bus0,sn]._body -=(r*network.model.passive_branch_p_sq_in[bt,bn,sn])
            network.model.power_balance[bus1,sn]._body -=(r*network.model.passive_branch_p_sq_out[bt,bn,sn])
    l_constraint(network.model,"mccormick_envelope_p", envelope,approx_index, list(passive_branches.index),snapshots,envelope_index)
    l_constraint(network.model,"losses", losses, list(passive_branches.index),snapshots)
    l_constraint(network.model,"envelope_intervalls",intervalls,approx_index, list(range(nrOfintervalls+1)), list(passive_branches.index),snapshots)
                                                