import pypsa
import numpy as np
import pandas as pd
import pyomo.environ as penv

from pypsa.opt import LConstraint, l_constraint, LExpression
from pyomo.environ import Var


def post_cosine(network, snapshots, duals):

    num_intervals = 10
    
    passive_branches = network.passive_branches()

    loss = pd.DataFrame(0, index=snapshots, columns=network.lines.index)

    for branch in passive_branches.index:
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]

        x = passive_branches.at[branch,"x_pu_eff"]
        b = 1/x

        if passive_branches.at[branch,"s_nom_extendable"]:
            xU = passive_branches.at[branch,"s_nom_max"] * x
        else:
            xU = passive_branches.at[branch,"s_nom"] * x
        
        for sn in snapshots:
            for i in range(num_intervals): 
                lower = xU * i / num_intervals
                max_val = lower + ( xU / num_intervals )
                slope = 2 * b * np.sin(lower)
                
                loss.loc[sn,bt] = 0.5 * slope * (
                                    network.model.delta_angle_positive[bt,bn,sn,i] + 
                                    network.model.delta_angle_negative[bt,bn,sn,i]
                                  )
                
                loss.loc[sn,bn] = 0.5 * slope * (
                                    network.model.delta_angle_positive[bt,bn,sn,i] +
                                    network.model.delta_angle_negative[bt,bn,sn,i]
                                  )    
    
    network.lines_t["loss"] = loss


# https://www.iit.comillas.edu/aramos/papers/losses.pdf
def cosine(network, snapshots):

    num_intervals = 10

    passive_branches = network.passive_branches()
    
    network.model.delta_angle_positive = Var(list(passive_branches.index),
                                             snapshots, list(range(num_intervals)))
    network.model.delta_angle_negative = Var(list(passive_branches.index),
                                             snapshots, list(range(num_intervals)))
    network.model.difference_constraint = penv.ConstraintList()
    
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
            xU = passive_branches.at[branch,"s_nom_max"] * x
        else:
            xU = passive_branches.at[branch,"s_nom"] * x

        for sn in snapshots:
            for i in range(num_intervals): 
                lower = xU * i / num_intervals
                max_val = lower + ( xU / num_intervals )
                slope = 2 * b * np.sin(lower)

                lhs = LExpression([(1,network.model.delta_angle_positive[bt,bn,sn,i])])
                intervals_lower[bt,bn,sn,i] = LConstraint(lhs, ">=", LExpression())

                lhs = LExpression([(1,network.model.delta_angle_negative[bt,bn,sn,i])],(-max_val))
                intervals_upper[bt,bn,sn,i] = LConstraint(lhs, "<=", LExpression())

                loss_term = 0.5 * slope * ( network.model.delta_angle_positive[bt,bn,sn,i] + \
                                            network.model.delta_angle_negative[bt,bn,sn,i] )
                network.model.power_balance[bus0,sn]._body -= loss_term
                network.model.power_balance[bus1,sn]._body -= loss_term
            
            network.model.difference_constraint.add(
                ( network.model.voltage_angles[bus0,sn] - \
                  network.model.voltage_angles[bus1,sn] ) == \
                ( sum( ( network.model.delta_angle_positive[bt,bn,sn,i] - \
                         network.model.delta_angle_negative[bt,bn,sn,i] ) \
                       for i in range(num_intervals))))
    
    l_constraint(network.model, "intervals_lower", intervals_lower,
                 list(passive_branches.index), snapshots, list(range(num_intervals)))

    l_constraint(network.model, "intervals_upper", intervals_upper,
                 list(passive_branches.index), snapshots, list(range(num_intervals)))   


def post_quadratic(network, snapshots, duals):
    
    passive_branches = network.passive_branches()

    loss = pd.DataFrame(0, index=snapshot, columns=network.lines.index)

    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch, "bus0"]
        bus1 = passive_branches.at[branch, "bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]
        
        x = passive_branches.at[branch,"x_pu_eff"]
        b = 1/x
        
        for sn in snapshots: 
            loss.loc[sn,bn] = (b*network.model.delta_angle_sq[bt,bn,sn])
    
    network.lines_t["loss"] =  loss

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
            xU = passive_branches.at[branch,"s_nom_max"] * x
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
            lhs = LExpression([(1,network.model.delta_angle_sq[bt,bn,sn])])
            rhs = LExpression([(xU,network.model.delta_angle[bt,bn,sn]),
                               (-xU,network.model.delta_angle[bt,bn,sn])],
                              (xU**2))
            envelope[bt,bn,sn,"wov1"] = LConstraint(lhs, "<=", rhs)

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


def post_lldc(network, snapshots, duals):

    passive_branches = network.passive_branches()
    
    loss = pd.DataFrame(0, index=snapshots, columns=network.lines.index)
    
    for branch in passive_branches.index:
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]
        
        r = passive_branches.at[branch, "r_pu_eff"]
        
        for sn in snapshots: 
            loss.loc[sn,bn] = r * network.model.passive_branch_p_sq_out[bt,bn,sn].value + \
                              r * network.model.passive_branch_p_sq_in[bt,bn,sn].value
    
    network.lines_t["loss"] = loss
        

#https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6345342
def lldc(network, snapshots):
    
    num_intervals = 10

    passive_branches = network.passive_branches()
    
    network.model.passive_branch_p_out = Var(list(passive_branches.index), snapshots, bounds=(0,float('inf')))
    network.model.passive_branch_p_in = Var(list(passive_branches.index), snapshots, bounds=(0,float('inf')))
    network.model.passive_branch_p_sq_in = Var(list(passive_branches.index), snapshots, bounds=(0,float('inf')))
    network.model.passive_branch_p_sq_out = Var(list(passive_branches.index), snapshots, bounds=(0,float('inf')))
    
    approx_index = ["in","out"]
    envelope_index = ["xup","wov1"]
    
    losses = {}
    envelope = {}
    intervals = {}
    
    for branch in passive_branches.index:
        bus0 = passive_branches.at[branch,"bus0"]
        bus1 = passive_branches.at[branch,"bus1"]
        sub = passive_branches.at[branch,"sub_network"]
        bt = branch[0]
        bn = branch[1]

        r = passive_branches.at[branch,"r_pu_eff"]

        if passive_branches.at[branch,"s_nom_extendable"]:
            xU = passive_branches.at[branch,"s_nom_max"]
        else:
            xU = passive_branches.at[branch,"s_nom"]
        
        xL = 0

        for sn in snapshots:
            lhs = LExpression([(1,network.model.passive_branch_p_out[bt,bn,sn]),
                               (-1,network.model.passive_branch_p_in[bt,bn,sn]),
                               (-1,network.model.passive_branch_p[bt,bn,sn])])
            losses[bt,bn,sn] = LConstraint(lhs, "==", LExpression())
            
            # bounds
            lhs = LExpression([(1,network.model.passive_branch_p_in[bt,bn,sn]),
                               (r,network.model.passive_branch_p_sq_in[bt,bn,sn]),
                               (r,network.model.passive_branch_p_sq_out[bt,bn,sn])],
                              (-xU))
            envelope["in",bt,bn,sn,"xup"] = LConstraint(lhs, "<=", LExpression())

            lhs = LExpression([(1,network.model.passive_branch_p_out[bt,bn,sn]),
                               (r,network.model.passive_branch_p_sq_in[bt,bn,sn]),
                               (r,network.model.passive_branch_p_sq_out[bt,bn,sn])],
                              (-xU))
            envelope["out",bt,bn,sn,"xup"] = LConstraint(lhs, "<=", LExpression())
            
            # approximation from above     
            lhs = LExpression([(1,network.model.passive_branch_p_sq_in[bt,bn,sn])])
            rhs = LExpression([(xU,network.model.passive_branch_p_in[bt,bn,sn]),
                               (xL,network.model.passive_branch_p_in[bt,bn,sn])],
                              (-xU*xL))
            envelope["in",bt,bn,sn,"wov1"] = LConstraint(lhs, "<=", rhs)

            lhs = LExpression([(1,network.model.passive_branch_p_sq_out[bt,bn,sn])])
            rhs = LExpression([(xU,network.model.passive_branch_p_out[bt,bn,sn]),
                               (xL,network.model.passive_branch_p_out[bt,bn,sn])],
                              (-xU*xL))
            envelope["out",bt,bn,sn,"wov1"] = LConstraint(lhs, "<=", rhs)
            
            # appoximation from below in 10 intervals
            for i in range(num_intervals+1):
                lower = xU * i / num_intervals

                lhs = LExpression([(1,network.model.passive_branch_p_sq_in[bt,bn,sn])])
                rhs = LExpression([(2*lower,network.model.passive_branch_p_in[bt,bn,sn])],(-lower**2))
                intervals["in",i,bt,bn,sn] = LConstraint(lhs, ">=", rhs)
                
                lhs = LExpression([(1,network.model.passive_branch_p_sq_out[bt,bn,sn])])
                rhs = LExpression([(2*lower,network.model.passive_branch_p_out[bt,bn,sn])],(-lower**2))
                intervals["out",i,bt,bn,sn] = LConstraint(lhs, ">=", rhs)
            
            # add p_sq to nodal power balance
            # use of ._body because of pyomo bug
            network.model.power_balance[bus0,sn]._body -= r * network.model.passive_branch_p_sq_in[bt,bn,sn]
            network.model.power_balance[bus1,sn]._body -= r * network.model.passive_branch_p_sq_out[bt,bn,sn]
   
    l_constraint(network.model, "mccormick_envelope_p", envelope,approx_index,
                 list(passive_branches.index), snapshots, envelope_index)

    l_constraint(network.model, "losses", losses,
                 list(passive_branches.index), snapshots)

    l_constraint(network.model, "envelope_intervals", intervals, approx_index,
                 list(range(num_intervals+1)), list(passive_branches.index), snapshots)
