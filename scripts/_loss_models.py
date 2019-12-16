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
                
    network.model.loss = Var(list(passive_branches.index), snapshots)
    
    intervals_upper = {}
    intervals_lower = {}
    
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch,"bus0"]
        bus1 = passive_branches.at[branch,"bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]
        
        x = passive_branches.at[branch,"x_pu_eff"]
        b = 1/x
        
        if passive_branches.at[branch,"s_nom_extendable"]:
            xU = passive_branches.at[branch,"s_nom_max"] * x * passive_branches.at[branch, "s_max_pu"]
        else:
            xU = passive_branches.at[branch,"s_nom"] * x 

        for sn in snapshots:
            for i in range(num_intervals): 
                lower = xU * i / num_intervals
                max_val = lower + ( xU / num_intervals )
                slope = 2 * b * np.sin(lower)

                lhs = LExpression([(1,network.model.delta_angle[bt,bn,sn,i])])
                intervals_lower[bt,bn,sn,i] = LConstraint(lhs, ">=", LExpression())

                lhs = LExpression([(1,network.model.delta_angle[bt,bn,sn,i])],(-max_val))
                intervals_upper[bt,bn,sn,i] = LConstraint(lhs, "<=", LExpression())

                network.model.power_balance[bus0,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
                network.model.power_balance[bus1,sn]._body -= 0.5 * network.model.loss[bt,bn,sn]
            
            lhs = LExpression([(1, network.model.delta_angle[bt,bn,sn])])
            rhs = LExpression([(1, network.model.voltage_angles[bus0,sn]),(-1, network.model.voltage_angles[bus1,sn])])
            difference[bt,bn,sn] = LConstraint(lhs, "=", rhs)

            lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
            rhs = LExpression([(slope,network.model.delta_angle[bt,bn,sn])])
            losses[bt,bn,sn] = LConstraint(lhs, "=", rhs)

    
    l_constraint(network.model, "intervals_lower", intervals_lower,
                 list(passive_branches.index), snapshots, list(range(num_intervals)))

    l_constraint(network.model, "intervals_upper", intervals_upper,
                 list(passive_branches.index), snapshots, list(range(num_intervals)))   

    l_constraint(network.model, "delta_angle", difference, list(passive_branches.index), snapshots)

    l_constraint(network.model, "losses", losses, list(passive_branches.index), snapshots)


# https://www.iit.comillas.edu/aramos/papers/losses.pdf
def quadratic(network, snapshots):
    
    num_intervals = 10

    passive_branches = network.passive_branches()
    
    network.model.delta_angle = Var(list(passive_branches.index), snapshots)
    network.model.delta_angle_sq = Var(list(passive_branches.index), snapshots)
    
    envelope_index = ["wov1", "xup", "xlow"]
    
    differences = {}
    envelope = {}
    intervals = {}
    
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
            xU = passive_branches.at[branch,"s_nom_max"] * x * passive_branches_at[branch,"s_max_pu"]
        else:
            xU = passive_branches.at[branch,"s_nom"] * x

        for sn in snapshots:
            lhs = LExpression([
                    (1,network.model.voltage_angles[bus0,sn]),
                    (-1,network.model.voltage_angles[bus1,sn]),
                    (-1,network.model.delta_angle[bt,bn,sn])
                  ])
            differences[bt,bn,sn] = LConstraint(lhs, "==", LExpression())
            
            # bounds 
            lhs = LExpression([(1,network.model.delta_angle[bt,bn,sn])],-xU)
            envelope[bt,bn,sn,"xup"] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1,network.model.delta_angle[bt,bn,sn])],xU)
            envelope[bt,bn,sn,"xlow"] = LConstraint(lhs, "<=", LExpression())                                                       

            # approximation from above
            lhs = LExpression([(1,network.model.delta_angle_sq[bt,bn,sn])], (-xU**2))
            envelope[bt,bn,sn,"wov1"] = LConstraint(lhs, "<=", LExpression())

            # appoximation from below in intervals
            for i in range(num_intervals+1):
                lower = xU * i / num_intervals

                lhs = LExpression([(1,network.model.delta_angle_sq[bt,bn,sn])])
                rhs = LExpression([(2*lower,network.model.delta_angle[bt,bn,sn])],(-lower**2))
                intervals[i,bt,bn,sn] = LConstraint(lhs, ">=" ,rhs)

            # add p_sq to nodal power balance
            # use of ._body because of pyomo bug
            network.model.power_balance[bus0,sn]._body -= 0.5 * g * network.model.delta_angle_sq[bt,bn,sn]
            network.model.power_balance[bus1,sn]._body -= 0.5 * g * network.model.delta_angle_sq[bt,bn,sn]
    
    l_constraint(network.model, "mccormick_envelope_delta_angle", envelope,
                 list(passive_branches.index),snapshots,envelope_index)

    l_constraint(network.model, "differences", differences,
                 list(passive_branches.index),snapshots)

    l_constraint(network.model, "envelope_intervals", intervals,
                 list(range(num_intervals+1)), list(passive_branches.index), snapshots)    
        

#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6345342
def lldc(network, snapshots):
    
    num_intervals = 10

    passive_branches = network.passive_branches()
    network.model.loss = Var(list(passive_branches.index), snapshots, bounds=(0,float('inf')))
    network.model.passive_branch_p_sq = Var(list(passive_branches.index), snapshots, bounds=(0,float('inf')))
    
    losses = {}
    upper_bound = {}
    intervals = {}
    
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch,"bus0"]
        bus1 = passive_branches.at[branch,"bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]

        r = passive_branches.at[branch,"r_pu_eff"]

        if passive_branches.at[branch,"s_nom_extendable"]:
            xU = passive_branches.at[branch,"s_nom_max"] * passive_branches_at * passive_branches.at[branch,"s_max"]
        else:
            xU = passive_branches.at[branch,"s_nom"]
        
        xL = 0

        for sn in snapshots:
            lhs = LExpression([(1,network.model.loss[bt,bn,sn])])
            rhs = LExpression([(r,network.model.passive_branch_p_sq[bt,bn,sn])])
            losses[bt,bn,sn] = LConstraint(lhs, "=", rhs)

            # approximation from above     
            lhs = LExpression([(1,network.model.passive_branch_p_sq[bt,bn,sn])])
            rhs = LExpression([(xU,network.model.passive_branch_p[bt,bn,sn]),
                               (xL,network.model.passive_branch_p[bt,bn,sn])],
                              (-xU*xL))
            upper_bound[bt,bn,sn] = LConstraint(lhs, "<=", rhs)

            
            # appoximation from below in 10 intervals
            for i in range(num_intervals+1):
                lower = xU * i / num_intervals

                lhs = LExpression([(1,network.model.passive_branch_p_sq[bt,bn,sn])])
                rhs = LExpression([(2*lower,network.model.passive_branch_p[bt,bn,sn])],(-lower**2))
                intervals[i,bt,bn,sn] = LConstraint(lhs, ">=", rhs)
                          
            # add p_sq to nodal power balance
            # use of ._body because of pyomo bug
            network.model.power_balance[bus0,sn]._body -= network.model.loss[bt,bn,sn]
            network.model.power_balance[bus1,sn]._body -= network.model.loss[bt,bn,sn]
   
    l_constraint(network.model, "upper_bound", upper_bound, list(passive_branches.index), snapshots)

    l_constraint(network.model, "losses", losses,
                 list(passive_branches.index), snapshots)

    l_constraint(network.model, "tangents", intervals, list(range(num_intervals+1)), list(passive_branches.index), snapshots)



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