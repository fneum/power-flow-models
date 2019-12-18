import pypsa
import numpy as np
import pandas as pd
import pyomo.environ as penv

from pypsa.opt import LConstraint, l_constraint, LExpression
from pyomo.environ import Var


# https://www.iit.comillas.edu/aramos/papers/losses.pdf
def cosine(network, snapshots):

    num_intervals = 10

    passive_branches = network.passive_branches()
    
    network.model.delta_angle = Var(list(passive_branches.index), snapshots)
                
    network.model.loss = Var(list(passive_branches.index), snapshots,)
    
    loss_upper = {}
    loss_lower = {}
    loss_tangents_neg = {}
    loss_tangents_pos = {}
    delta_lower = {}
    delta_upper = {}
    difference = {}

    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch,"bus0"]
        bus1 = passive_branches.at[branch,"bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]
        
        x = passive_branches.at[branch,"x_pu_eff"]
        b = 1/x
        r = passive_branches.at[branch,"r_pu_eff"]
        g = 1/r
        
        if passive_branches.at[branch,"s_nom_extendable"]:
            xU = passive_branches.at[branch,"s_nom_max"] * x * passive_branches.at[branch, "s_max_pu"]
        else:
            xU = passive_branches.at[branch,"s_nom"] * x  * passive_branches.at[branch, "s_max_pu"]
        max_loss = 2 * g * (1 - np.cos(xU))
        xL = 0
        for sn in snapshots:
            for i in range(num_intervals +1): 
                lower = xU * (i / num_intervals)
                max_val = lower + ( xU / num_intervals )
                slope =  g * np.sin(lower)
                t = (2 * g) * (1 - np.cos(lower)) - (slope * lower)

                lhs = LExpression([(1, network.model.loss[bt,bn,sn])])
                rhs = LExpression([(slope,network.model.delta_angle[bt,bn,sn])], (-t))
                loss_tangents_pos[bt,bn,sn,i] = LConstraint(lhs, ">=", rhs)

                lhs = LExpression([(1, network.model.loss[bt,bn,sn])])
                rhs = LExpression([(-slope, network.model.delta_angle[bt,bn,sn])], (-t))
                loss_tangents_neg[bt,bn,sn,i] = LConstraint(lhs, ">=", rhs)

            network.model.power_balance[bus0,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.power_balance[bus1,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.flow_lower[bt,bn,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.flow_upper[bt,bn,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]

            lhs = LExpression([(1,network.model.delta_angle[bt,bn,sn])],-xU)
            delta_lower[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1,network.model.delta_angle[bt,bn,sn])],xU)
            delta_upper[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())  
            
            lhs = LExpression([(1, network.model.loss[bt,bn,sn])],(-max_loss))
            loss_upper[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
            loss_lower[bt,bn,sn] = LConstraint(lhs, ">=", LExpression())

            lhs = LExpression([(1, network.model.delta_angle[bt,bn,sn])])
            rhs = LExpression([(1, network.model.voltage_angles[bus0,sn]),(-1, network.model.voltage_angles[bus1,sn])])
            difference[bt,bn,sn] = LConstraint(lhs, "==", rhs)

    
    l_constraint(network.model, "loss_tangents_neg", loss_tangents_neg, list(passive_branches.index), snapshots, list(range(num_intervals + 1)))

    l_constraint(network.model, "loss_tangents_pos", loss_tangents_pos, list(passive_branches.index), snapshots, list(range(num_intervals + 1)))

    l_constraint(network.model, "loss_upper", loss_upper, list(passive_branches.index), snapshots)

    l_constraint(network.model, "loss_lower", loss_lower, list(passive_branches.index), snapshots)   

    l_constraint(network.model, "angle_difference", difference, list(passive_branches.index), snapshots)

    l_constraint(network.model, "delta_upper", delta_upper, list(passive_branches.index), snapshots)

    l_constraint(network.model, "delta_lower", delta_lower, list(passive_branches.index), snapshots)



# https://www.iit.comillas.edu/aramos/papers/losses.pdf
def quadratic(network, snapshots):
    
    num_intervals = 10

    passive_branches = network.passive_branches()
    
    network.model.delta_angle = Var(list(passive_branches.index), snapshots)
    network.model.loss = Var(list(passive_branches.index), snapshots)
        
    differences = {}
    loss_lower = {}
    loss_upper = {}
    loss_tangents_pos = {}
    loss_tangents_neg = {}
    delta_lower = {}
    delta_upper = {}
    
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch,"bus0"]
        bus1 = passive_branches.at[branch,"bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]

        x = passive_branches.at[branch,"x_pu_eff"]
        b = 1/x
        r = passive_branches.at[branch,"r_pu_eff"]
        g = 1/r
        
        if passive_branches.at[branch,"s_nom_extendable"]:
            xU = passive_branches.at[branch,"s_nom_max"] * x * passive_branches.at[branch,"s_max_pu"]
        else:
            xU = passive_branches.at[branch,"s_nom"] * x * passive_branches.at[branch, "s_max_pu"]

        for sn in snapshots:
            lhs = LExpression([
                    (1,network.model.voltage_angles[bus0,sn]),
                    (-1,network.model.voltage_angles[bus1,sn]),
                    (-1,network.model.delta_angle[bt,bn,sn])
                  ])
            differences[bt,bn,sn] = LConstraint(lhs, "==", LExpression())
            
            # bounds 
            lhs = LExpression([(1,network.model.delta_angle[bt,bn,sn])],-xU)
            delta_lower[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1,network.model.delta_angle[bt,bn,sn])],xU)
            delta_upper[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())  

            lhs = LExpression([(1,network.model.loss[bt,bn,sn])], -g*(xU**2))
            loss_upper[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
            loss_lower[bt,bn,sn] = LConstraint(lhs, ">=", LExpression())                                                     

            # approximation from above  :
            #is this necessary would replace loss_upper ??
            #lhs = LExpression([(1,network.model.loss[bt,bn,sn])], (-r*(xU**2)))
            #envelope[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            # appoximation from below in intervals
            for i in range(num_intervals+1):
                lower = xU * i / num_intervals

                lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
                rhs = LExpression([(2*g*lower,network.model.delta_angle[bt,bn,sn])],(-g*lower**2))
                loss_tangents_pos[bt,bn,sn, i] = LConstraint(lhs, ">=" ,rhs)

                lhs = LExpression([(1, network.model.loss[bt,bn,sn])])
                rhs = LExpression([(2*g*lower, network.model.delta_angle[bt,bn,sn])], (-g*lower**2))
                loss_tangents_neg[bt,bn,sn, i] = LConstraint(lhs, ">=", rhs)

            # add p_sq to nodal power balance
            # use of ._body because of pyomo bug
            network.model.power_balance[bus0,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.power_balance[bus1,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.flow_lower[bt,bn,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.flow_upper[bt,bn,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
    
#    l_constraint(network.model, "mccormick_envelope_delta_angle", envelope,
 #                list(passive_branches.index),snapshots,envelope_index)

    l_constraint(network.model, "loss_tangents_neg", loss_tangents_neg, list(passive_branches.index), snapshots, list(range(num_intervals + 1)))

    l_constraint(network.model, "loss_tangents_pos", loss_tangents_pos, list(passive_branches.index), snapshots, list(range(num_intervals + 1)))

    l_constraint(network.model, "loss_upper", loss_upper, list(passive_branches.index), snapshots)

    l_constraint(network.model, "loss_lower", loss_lower, list(passive_branches.index), snapshots)   

    l_constraint(network.model, "angle_difference", differences, list(passive_branches.index), snapshots)

    l_constraint(network.model, "delta_upper", delta_upper, list(passive_branches.index), snapshots)

    l_constraint(network.model, "delta_lower", delta_lower, list(passive_branches.index), snapshots)
        

#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6345342
def lldc(network, snapshots):
    
    num_intervals = 10

    passive_branches = network.passive_branches()
    network.model.loss = Var(list(passive_branches.index), snapshots)
    
    upper_bound = {}
    loss_upper = {}
    loss_lower = {}
    loss_tangents_neg = {}
    loss_tangents_pos = {}
    
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch,"bus0"]
        bus1 = passive_branches.at[branch,"bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]

        r = passive_branches.at[branch,"r_pu_eff"]

        if passive_branches.at[branch,"s_nom_extendable"]:
            xU = passive_branches.at[branch,"s_nom_max"] * x * passive_branches.at[branch,"s_max"]
        else:
            xU = passive_branches.at[branch,"s_nom"] * passive_branches.at[branch, "s_max_pu"]
        
        xL = -xU

        for sn in snapshots:
            # approximation from above     
            #lhs = LExpression([(r,network.model.loss[bt,bn,sn])])
            #rhs = LExpression([(r*xU,network.model.passive_branch_p[bt,bn,sn]),
            #                   (r*xL,network.model.passive_branch_p[bt,bn,sn])],
            #                  (-r*xU*xL))
            #upper_bound[bt,bn,sn] = LConstraint(lhs, "<=", rhs)

            lhs = LExpression([(1,network.model.loss[bt,bn,sn])], -r*(xU**2))
            loss_upper[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
            loss_lower[bt,bn,sn] = LConstraint(lhs, ">=", LExpression())      

            # appoximation from below in 10 intervals
            for i in range(num_intervals+1):
                lower = xU * i / num_intervals

                lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
                rhs = LExpression([((r*2*lower),network.model.passive_branch_p[bt,bn,sn])],(-r*(lower**2)))
                loss_tangents_pos[i,bt,bn,sn] = LConstraint(lhs, ">=", rhs)

                lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
                rhs = LExpression([((-r*2*lower), network.model.passive_branch_p[bt,bn,sn])],(-r*(lower**2)))
                loss_tangents_neg[i,bt,bn,sn] = LConstraint(lhs, ">=", rhs)
                          
            # add p_sq to nodal power balance
            # use of ._body because of pyomo bug
            network.model.power_balance[bus0,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.power_balance[bus1,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.flow_lower[bt,bn,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.flow_upper[bt,bn,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
   
    #l_constraint(network.model, "upper_bound", upper_bound, list(passive_branches.index), snapshots)

    l_constraint(network.model, "loss_upper", loss_upper, list(passive_branches.index), snapshots)

    l_constraint(network.model, "loss_lower", loss_lower, list(passive_branches.index), snapshots)  

    l_constraint(network.model, "loss_tangents_neg", loss_tangents_pos, list(range(num_intervals+1)), list(passive_branches.index), snapshots)

    l_constraint(network.model, "loss_tangents_pos", loss_tangents_neg, list(range(num_intervals+1)), list(passive_branches.index),snapshots)


def post_processing(network, snapshots, duals):

    passive_branches = network.passive_branches()
    
    loss = pd.DataFrame(0, index=snapshots, columns=network.lines.index)
    
    for branch in passive_branches.index:
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]
        
        r = passive_branches.at[branch, "r_pu_eff"]
        
        for sn in snapshots: 

            loss.loc[sn,bt] = 0.5 * network.model.loss[bt,bn,sn]
            loss.loc[sn,bn] = 0.5 * network.model.loss[bt,bn,sn]
    network.lines_t["loss"] = loss