import pypsa
import numpy as np
import pandas as pd

from pypsa.opt import LConstraint, l_constraint, LExpression
from pyomo.environ import Var, value


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
        bt = branch[0]
        bn = branch[1]
        
        x_pu_eff = passive_branches.at[branch,"x_pu_eff"]
        x = passive_branches.at[branch,"x"]
        r_pu_eff = passive_branches.at[branch,"r_pu_eff"]
        r = passive_branches.at[branch,"r"]
        g = r / (x**2)
        
        if passive_branches.at[branch,"s_nom_extendable"]:
            d_theta_max = passive_branches.at[branch,"s_nom_max"] * x_pu_eff * passive_branches.at[branch, "s_max_pu"]
        else:
            d_theta_max = passive_branches.at[branch,"s_nom"] * x_pu_eff  * passive_branches.at[branch, "s_max_pu"]
        
        max_loss = 2 * g * (1 - np.cos(d_theta_max))
        
        for sn in snapshots:
            for i in range(num_intervals+1):
                d_theta_i = d_theta_max * (i / num_intervals)
                losses_i = 2 * g * (1 - np.cos(d_theta_i)) 
                slope_i =  2 * g * np.sin(d_theta_i)
                offset_i = losses_i - (slope_i * d_theta_i)

                lhs = LExpression([(1, network.model.loss[bt,bn,sn])])
                rhs = LExpression([(slope_i,network.model.delta_angle[bt,bn,sn])], offset_i)
                loss_tangents_pos[bt,bn,sn,i] = LConstraint(lhs, ">=", rhs)

                lhs = LExpression([(1, network.model.loss[bt,bn,sn])])
                rhs = LExpression([(-slope_i, network.model.delta_angle[bt,bn,sn])], offset_i)
                loss_tangents_neg[bt,bn,sn,i] = LConstraint(lhs, ">=", rhs)

            network.model.power_balance[bus0,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.power_balance[bus1,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.flow_lower[bt,bn,sn]._body -= network.model.loss[bt,bn,sn]
            network.model.flow_upper[bt,bn,sn]._body -= network.model.loss[bt,bn,sn]

            lhs = LExpression([(1, network.model.delta_angle[bt,bn,sn])], (-d_theta_max))
            delta_upper[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1, network.model.delta_angle[bt,bn,sn])], d_theta_max)
            delta_lower[bt,bn,sn] = LConstraint(lhs, ">=", LExpression())  
            
            lhs = LExpression([(1, network.model.loss[bt,bn,sn])], (-max_loss))
            loss_upper[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1, network.model.loss[bt,bn,sn])])
            loss_lower[bt,bn,sn] = LConstraint(lhs, ">=", LExpression())

            lhs = LExpression([(1, network.model.delta_angle[bt,bn,sn])])
            rhs = LExpression([(1, network.model.voltage_angles[bus0,sn]),
                               (-1, network.model.voltage_angles[bus1,sn])])
            difference[bt,bn,sn] = LConstraint(lhs, "==", rhs)

    
    l_constraint(network.model, "loss_tangents_neg", loss_tangents_neg,
                 list(passive_branches.index), snapshots, list(range(num_intervals + 1)))

    l_constraint(network.model, "loss_tangents_pos", loss_tangents_pos,
                 list(passive_branches.index), snapshots, list(range(num_intervals + 1)))

    l_constraint(network.model, "loss_upper", loss_upper, list(passive_branches.index), snapshots)

    l_constraint(network.model, "loss_lower", loss_lower, list(passive_branches.index), snapshots)   

    l_constraint(network.model, "angle_difference", difference, list(passive_branches.index), snapshots)

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
        bt = branch[0]
        bn = branch[1]

        r_pu_eff = passive_branches.at[branch,"r_pu_eff"]

        if passive_branches.at[branch,"s_nom_extendable"]:
            p_max = passive_branches.at[branch,"s_nom_max"] * passive_branches.at[branch,"s_max_pu"]
        else:
            p_max = passive_branches.at[branch,"s_nom"] * passive_branches.at[branch, "s_max_pu"]
        
        for sn in snapshots:

            lhs = LExpression([(1,network.model.loss[bt,bn,sn])], -r_pu_eff*(p_max**2))
            loss_upper[bt,bn,sn] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
            loss_lower[bt,bn,sn] = LConstraint(lhs, ">=", LExpression())      

            for i in range(num_intervals+1):
                p_i = p_max * i / num_intervals
                losses_i = r_pu_eff * p_i**2
                slope_i = r_pu_eff * 2 * p_i
                offset_i = losses_i - slope_i * p_i

                lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
                rhs = LExpression([(slope_i, network.model.passive_branch_p[bt,bn,sn])], offset_i)
                loss_tangents_pos[i,bt,bn,sn] = LConstraint(lhs, ">=", rhs)

                lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
                rhs = LExpression([(-slope_i, network.model.passive_branch_p[bt,bn,sn])], offset_i)
                loss_tangents_neg[i,bt,bn,sn] = LConstraint(lhs, ">=", rhs)
                          
            # add p_sq to nodal power balance
            # use of ._body because of pyomo bug
            network.model.power_balance[bus0,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.power_balance[bus1,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            network.model.flow_lower[bt,bn,sn]._body -= network.model.loss[bt,bn,sn]
            network.model.flow_upper[bt,bn,sn]._body -= network.model.loss[bt,bn,sn]
   
    l_constraint(network.model, "loss_upper", loss_upper, list(passive_branches.index), snapshots)

    l_constraint(network.model, "loss_lower", loss_lower, list(passive_branches.index), snapshots)  

    l_constraint(network.model, "loss_tangents_neg", loss_tangents_pos, list(range(num_intervals+1)), list(passive_branches.index), snapshots)

    l_constraint(network.model, "loss_tangents_pos", loss_tangents_neg, list(range(num_intervals+1)), list(passive_branches.index),snapshots)


def post_processing(network, snapshots, duals):

    passive_branches = network.passive_branches()
    
    loss = pd.DataFrame(0, index=snapshots, columns=network.lines.index)
    
    for branch in passive_branches.index:
        bt = branch[0]
        bn = branch[1]
        for sn in snapshots: 
            loss.loc[sn,bn] = value(network.model.loss[bt,bn,sn])
    
    network.lines_t["loss"] = loss