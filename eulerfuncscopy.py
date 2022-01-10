# this file implements the various step functions to propagate the relevant variables
# as well as the euler function itself and the graphing function

# TODO implement cleanup after function ends early i.e singularity, ensure outside points are set to 0

import json
import numpy as np
from sympy import *
from mayavi import mlab
from consts import *
import math
import pickle

u, v, alpha, beta, phi, theta, t, s = symbols(
    'u v alpha beta phi theta t s')


class Vec2D:
    # will represent the vector n or r at each point

    def __init__(self, x, y):
        self.x = x
        self.y = y

    # perp yields the perpendicular of this vector
    def perp(self):
        return Vec2D(-self.y, self.x)

    # normalize returns this vector normalized
    def normalize(self):
        mag = np.sqrt(float(self.x*self.x + self.y*self.y))
        if mag != 0:
            return Vec2D(self.x/mag, self.y/mag)
        else:
            return self

    def add(self, vec):
        return Vec2D(self.x + vec.x, self.y + vec.y)

    def times(self, num):
        return Vec2D(self.x*num, self.y*num)


class Vec3D:
    # will represent the vector ra

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    # normalize returns this vector normalized
    def normalize(self):
        mag = np.sqrt(float(self.x**2 + self.y**2 + self.z**2))
        if mag != 0:
            return Vec3D(self.x/mag, self.y/mag, self.z/mag)
        else:
            return self

    def add(self, vec):
        return Vec3D(self.x + vec.x, self.y + vec.y, self.z + vec.z)

    def times(self, num):
        return Vec3D(self.x*num, self.y*num, self.z*num)


class DataPoint:
    # point object, contains all necesary values at each point,
    # note uv represents distance from U0,V0
    def __init__(self, uv, b, s, beta, alpha, rs, ns, Ks, epsilon, p, dup, q):
        self.uv = uv
        self.b = b
        self.s = s
        self.beta = beta
        self.alpha = alpha
        # r and n 2D vectors
        self.r = rs
        self.n = ns
        self.K = Ks  # curvature on final surface at point
        self.epsilon = epsilon
        self.p = p
        # u deriv of p:
        self.dup = dup
        self.q = q


def list_dromegadv(oldpoint, Lambda):
    # Lambda = vector of lambdas e.g. [0,1,2]
    out = []
    for i in range(len(oldpoint.r)):
        betaa = oldpoint.beta*LAMBDA2(Lambda[i], oldpoint)
        nomega = oldpoint.n[i]
        out.append(Vec2D(betaa*(nomega.perp().x), betaa*(nomega.perp().y)))
    return out


def list_dnomegadv(mathobj, Data, uold, vold, Lambda):
    # perhaps fix so uses new and old point rather than old and older point?
    # so maybe use oldpoint for beta etc, and newpoint for epsilon when doing manual difference
    oldpoint = Data[vold][uold]
    out = []
    #print(uold, " ", vold)
    for i in range(len(oldpoint.n)):
        mult = oldpoint.beta * LAMBDA2(Lambda[i], oldpoint)
        dlambda_term = (LAMBDA2(Lambda[i], oldpoint) - LAMBDA2(
            Lambda[i], Data[vold - 1][uold]))/(oldpoint.epsilon - Data[vold - 1][uold].epsilon)
        dlambda_term *= (oldpoint.p /
                         (LAMBDA2(Lambda[i], oldpoint)*LAMBDA1(Lambda[i], oldpoint)))
        kgv = oldpoint.s/LAMBDA1(Lambda[i], oldpoint) + dlambda_term
        first_term = Vec2D(oldpoint.n[i].perp().x,
                           oldpoint.n[i].perp().y).times(kgv)

        nprod = Matrix([[(oldpoint.n[i].perp()).x * oldpoint.n[i].x, (oldpoint.n[i].perp()).x * oldpoint.n[i].y],
                        [(oldpoint.n[i].perp()).y * oldpoint.n[i].x, (oldpoint.n[i].perp()).y * oldpoint.n[i].y]])
        G_u = mathobj.christoffel(oldpoint, i, 1)
        G_v = mathobj.christoffel(oldpoint, i, 2)
        u_term = G_u[0, 0]*nprod[0, 0] + G_u[0, 1]*nprod[0, 1] + \
            G_u[1, 0]*nprod[1, 0]+G_u[1, 1]*nprod[1, 1]
        v_term = G_v[0, 0]*nprod[0, 0] + G_v[0, 1]*nprod[0, 1] + \
            G_v[1, 0]*nprod[1, 0]+G_v[1, 1]*nprod[1, 1]
        second_term = (Vec2D(u_term, v_term))
        out.append(first_term.add(second_term.times(-1)).times(mult))
    return out


def dsdu_double(oldpoint, d):
    return oldpoint.alpha * d[1]


def dbetadu(oldpoint):
    s = oldpoint.s
    alpha = oldpoint.alpha
    beta = oldpoint.beta
    return s*alpha*beta


def dbdv_double(oldpoint, d):
    return oldpoint.beta * d[0]


def dalphadv(oldpoint):
    # same for multiple surfaces as before
    alpha = oldpoint.alpha
    beta = oldpoint.beta
    b = oldpoint.b
    return -alpha*beta*b


def dqdv(oldpoint, d):
    return oldpoint.beta*d[2]


def list_n_perp_step_double(mathobj, Data, u, v, sgn, d, Lambda):
    # now includes actuated!
    oldpoint = Data[v - 1][u]
    Data[v][u].b = oldpoint.b + sgn*STEP*dbdv_double(oldpoint, d)
    Data[v][u].alpha = oldpoint.alpha + sgn*STEP*dalphadv(oldpoint)
    deltar = list_dromegadv(oldpoint, Lambda)
    # deltan = list_dnomegadv(mathobj, Data, u, v-1, Lambda)
    # try from current point since already have epsilon - the only needed data at a point
    deltan = list_dnomegadv(mathobj, Data, u, v-1, Lambda)
    if u != 0:
        for i in range(len(oldpoint.n)):
            Data[v][u].r[i] = oldpoint.r[i].add(deltar[i].times(sgn*STEP))
            # print(Data[v][u].r[i].x, Data[v][u].r[i].y , " uv:", u, v, "\n")
            Data[v][u].n[i] = oldpoint.n[i].add(
                deltan[i].times(sgn*STEP)).normalize()
    for i in range(len(oldpoint.n)):
        rs_u = Data[v][u].r[i].x
        rs_v = Data[v][u].r[i].y
        Data[v][u].K[i] = mathobj.KGaussian(rs_u, rs_v, i)
    Data[v][u].q = oldpoint.q + sgn*STEP*dqdv(oldpoint, d)
    # forgot!! have to update p...
    # do it numerically?
    # p = du_e/alpha so do p=epsilon-old_u epsilon/ old alpha (or new alpha?)
    Data[v][u].p = (Data[v][u].epsilon - oldpoint.epsilon) / \
        (STEP*oldpoint.alpha)
    # seems loose to take double derivative numerically...
    Data[v][u].dup = (Data[v][u].p - oldpoint.p) / STEP


def n_step_double(Data, u, v, sgn, d):
    # still works for LIST
    oldpoint = Data[v][u-1]
    Data[v][u].beta = float(oldpoint.beta + sgn*STEP*dbetadu(oldpoint))
    Data[v][u].s = oldpoint.s + sgn*STEP*dsdu_double(oldpoint, d)


def at_perp_singularity(Data, u, v, sgn):
    # check if alpha is essentially 0, this hit a singularity
    # alpha updated on perp steps, so only check at last v
    # no longer check bounds since different for different surfaces
    # for i in range(len(Data[v][u].r)):
    #     if Data[v-1][u].r[i].x > BOUNDS or Data[v-1][u].r[i].y > BOUNDS:
    #         print("out of bounds")
    #         return True
    if Data[v-1][u].alpha < DELTA:
        print("alpha negative")
    return (Data[v-1][u].alpha < DELTA)


def at_parallel_singularity(Data, u, v, sgn):
    # check if beta is essentially 0, this hit a singularity
    # beta updated on parallel steps, so only check at last u
    # LIST check all rs
    # for i in range(len(Data[v][u].r)):
    #     if  Data[v][u].r[i].x > BOUNDS or Data[v][u-1].r[i].y > BOUNDS:
    #         print("out of bounds")
    #         return True
    if Data[v][u-1].beta < DELTA:
        print("beta negative")
    return (Data[v][u-1].beta < DELTA)


def generate_init_uv_numeric2(mathobj, Data, Lu, Lv, usgn, vsgn):
    # wholly numeric approach
    # init t and s
    su = 0
    tu = 0
    for i in range(1, THETASIZE):
        df = np.sqrt(float(Lu[0].subs(t, STEP)**2 + Lu[1].subs(t, STEP)**2))
        dt = STEP/(df)  # now stepping const ds
        su += STEP
        tu += dt
        deltax = (Lu[0].subs(t, tu) - Lu[0].subs(t, tu-dt))/dt
        deltay = (Lu[1].subs(t, tu) - Lu[1].subs(t, tu-dt))/dt
        oldpoint = Data[0][i-1]
        Data[0][i].uv = oldpoint.uv.add(
            Vec2D(deltax*usgn*STEP, deltay*usgn*STEP))
        for i in range(len(oldpoint.r)):
            Data[0][i].r[i] = oldpoint.r[i].add(
                Vec2D(deltax*LAMBDA1(Lambda[i], oldpoint)*usgn*STEP, deltay*LAMBDA1(Lambda[i], oldpoint)*usgn*STEP))
            Data[0][i].n[i] = Vec2D(deltax, deltay).normalize()
            Data[0][i].K[i] = mathobj.KGaussian(
                Data[0][i].r[i].x, Data[0][i].r[i].y, i)

        kgu = N(mathobj.geodesic_curvature(Data[0][i], 1, tu))
        Data[0][i].b = 1 * kgu  # maybe right?

    sv = 0
    tv = 0
    # need add list stuff
    # do same below:
    for i in range(1, PHISIZE):
        df = np.sqrt(float(Lv[0].subs(t, STEP)**2 + Lv[1].subs(t, STEP)**2))
        dt = STEP/(df)  # now stepping const ds
        sv += STEP
        tv += dt
        deltax = (Lv[0].subs(t, tv) - Lv[0].subs(t, tv-dt))/dt
        deltay = (Lv[1].subs(t, tv) - Lv[1].subs(t, tv-dt))/dt

        oldpoint = Data[i-1][0]
        Data[i][0].uv = oldpoint.uv.add(
            Vec2D(deltax*vsgn*STEP, deltay*vsgn*STEP))
        for i in range(len(oldpoint.r)):
            Data[i][0].r[i] = oldpoint.r[i].add(
                Vec2D(deltax*LAMBDA2(Lambda[i], oldpoint)*vsgn*STEP, deltay*LAMBDA2(Lambda[i], oldpoint)*vsgn*STEP))
            Data[i][0].n[i] = Vec2D(deltax, deltay).normalize()
            Data[i][0].K[i] = mathobj.KGaussian(
                Data[i][0].r[i].x, Data[i][0].r[i].y, i)

        kgv = N(mathobj.geodesic_curvature(Data[i][0], 2, tv))
        Data[i][0].s = 1 * kgv

    print("uv setup done")


def hit_nan(point):
    isnan = False
    # NOTE added checking if lambda at pt is less than 0...
    # is this a reasonable assumption that lambda > 0?
    # for i in range(3):
    #     if LAMBDA1(i, point) < 0 :
    #         print("Lamda1 ", i, " less than 0")
    #         isnan = True
    #     if LAMBDA2(i, point) < 0:
    #         print("Lamda2 ", i, " less than 0")
    #         isnan = True
    if math.isnan(point.b):
        print("b is nan")
        isnan = True
    if math.isnan(point.s):
        print("s is nan")
        isnan = True
    if math.isnan(point.beta):
        print("beta is nan")
        isnan = True
    if math.isnan(point.alpha):
        print("alpha is nan")
        isnan = True
    for i in range(len(point.K)):
        if math.isnan(point.K[i]):
            print("K ", i, " is nan")
            isnan = True
    if math.isnan(point.p):
        print("p is nan")
        isnan = True
    if math.isnan(point.epsilon):
        print("epsilon is nan")
        isnan = True
    if math.isnan(point.dup):
        print("dup is nan")
        isnan = True
    if math.isnan(point.q):
        print("q is nan")
        isnan = True
    if isnan:
        print("at v u: ", point.uv.y, " ", point.uv.x)
    return isnan


def gen_empty_Data(quadrant):
    # generates empty list of ~zeroed-out points
    u_sgn = 1 if (quadrant == 1 or quadrant == 4) else -1
    v_sgn = 1 if (quadrant == 1 or quadrant == 2) else -1
    return [[DataPoint(Vec2D(u*u_sgn*STEP, v*v_sgn*STEP), 0.0, 0.0, 0.0, 0.0, [Vec2D(0.0, 0.0), Vec2D(0, 0), Vec2D(0, 0)],
                       [Vec2D(0.0, 0.0), Vec2D(0, 0), Vec2D(0, 0)], [1.0, 1.0, 1.0], 1.0, 0.0, 0.0, 0.0) for u in range(THETASIZE)]
            for v in range(PHISIZE)]


def euler_procedure(mathobj, quadrant):
    # TODO stops halfway in u? e.g. if double bounds, u singularity doubles...
    # TODO plots, but not working with numeric uv init.
    u_sgn = 1 if (quadrant == 1 or quadrant == 4) else -1
    v_sgn = 1 if (quadrant == 1 or quadrant == 2) else -1

    # set up zeroed out data table and u,v lines:
    Data = [[DataPoint(Vec2D(u*u_sgn*STEP, v*v_sgn*STEP), 0.0, 0.0, 0.0, 0.0, [Vec2D(0.0, 0.0), Vec2D(0, 0), Vec2D(0, 0)],
                       [Vec2D(0.0, 0.0), Vec2D(0, 0), Vec2D(0, 0)], [1.0, 1.0, 1.0], 0.5, 0.0, 0.0, 0.0) for u in range(THETASIZE)]
            for v in range(PHISIZE)]
    # note, setting init epsilon = 0.5 above so no division by 0
    # i.e. ^ half material half cavity

    # init u and v curves:
    for i in range(1, THETASIZE):
        # u line
        u_pos = i*u_sgn*STEP*LAMBDA1(Lambda[0], Data[0][i-1])
        u1_3d_pos = i*u_sgn*STEP*LAMBDA1(Lambda[1], Data[0][i-1])
        u2_3d_pos = i*u_sgn*STEP*LAMBDA1(Lambda[2], Data[0][i-1])
        # good choice of epsilon, p, q?
        Data[0][i] = DataPoint(Vec2D(u_pos, 0.0), 0, 0, 0, 1,
                               [Vec2D(u_pos + U0, V0),
                                Vec2D(u1_3d_pos + U0, V0),
                                Vec2D(u2_3d_pos + U0, V0)],
                               [Vec2D(1.0, 0), Vec2D(1.0, 0), Vec2D(1.0, 0)],
                               [mathobj.KGaussian(u_pos + U0, V0, 0),
                                mathobj.KGaussian(u1_3d_pos + U0, V0, 1),
                                mathobj.KGaussian(u2_3d_pos + U0, V0, 2)],
                               0.5 - i*0.01, -0.01, 0.0, -0.01)  # q=0.0 doesn't seem to work... immediately cannot divide by 0
    # for any changing of signs i.e. u decreasing rather than increasing, eventually flips
    # so "wants" to increase along a line

    # cannot do this since no epsilon data, must calculate anew at each row in the main loop
    for i in range(1, PHISIZE):
        # v line
        v_pos = i*v_sgn*LAMBDA2(Lambda[0], Data[i-1][0])*STEP
        v1_3d_pos = i*v_sgn*LAMBDA2(Lambda[1], Data[i-1][0])*STEP
        v2_3d_pos = i*v_sgn*LAMBDA2(Lambda[2], Data[i-1][0])*STEP

        Data[i][0] = DataPoint(Vec2D(0.0, v_pos), 0, 0, 1, 0,
                               [Vec2D(U0, v_pos + V0),
                                Vec2D(U0, v1_3d_pos + V0),
                                Vec2D(U0, v2_3d_pos + V0)],
                               [Vec2D(1.0, 0), Vec2D(1.0, 0), Vec2D(1.0, 0)],
                               [mathobj.KGaussian(U0, v_pos + V0, 0),
                                mathobj.KGaussian(U0, v1_3d_pos + V0, 1),
                                mathobj.KGaussian(U0, v2_3d_pos + V0, 2)],
                               0.0, 0.0, 0.0, 0.0)

    # origin point:
    Data[0][0] = DataPoint(Vec2D(0.0, 0.0), 0, 0, 1, 1,
                           [Vec2D(U0, V0), Vec2D(U0, V0), Vec2D(U0, V0)],
                           [Vec2D(1.0, 0), Vec2D(1.0, 0), Vec2D(1.0, 0)],
                           [mathobj.KGaussian(U0, V0, 0), mathobj.KGaussian(
                               U0, V0, 1), mathobj.KGaussian(U0, V0, 2)],
                           0.5, -0.01, 0.0, -0.01)

    # generate_init_uv_numeric2(mathobj, Data, Matrix(
    #     [mathobj.curve1u, mathobj.curve1v]), Matrix([mathobj.curve2u, mathobj.curve2v]), u_sgn, v_sgn)

    # no longer need set b s manually when do init uv subroutine ^
    # integrate beta, s along u:
    for i in range(1, THETASIZE):
        oldpoint = Data[0][i-1]
        try:
            parallel_d = mathobj.d(Data[0][i-1], Data[0][i])
        except:
            print(
                "d not able to be calculated, likely M not invertible, i.e. hit singularity")
            return
        #parallel_d = mathobj.d(Data[0][i-1], Data[0][i])
        Data[0][i].s = oldpoint.s + u_sgn*STEP * \
            dsdu_double(oldpoint, parallel_d)
        Data[0][i].beta = oldpoint.beta + u_sgn*STEP*dbetadu(oldpoint)

    # see if any singularity along u inital line:
    u_singularity = 0
    for i in range(1, THETASIZE):
        if at_parallel_singularity(Data, i, 0, u_sgn):
            # record where singularity located:
            u_singularity = i
            break
    if u_singularity != 0:
        for i in range(u_singularity, THETASIZE):
            # zero out all points past the singularity:
            Data[0][i] = DataPoint(Vec2D(0, 0), 0.0, 0.0, 0.0, 0.0,
                                   [Vec2D(0.0, 0.0), Vec2D(0, 0), Vec2D(0, 0)],
                                   [Vec2D(0.0, 0.0), Vec2D(0, 0), Vec2D(0, 0)],
                                   [1.0, 1.0, 1.0], 0.0, 0.0, 0.0, 0.0)

    # now do iterative scheme:
    for v in range(1, PHISIZE):
        # update init v curve: THIS BREAKS IT BUT WHY? whoops did [i][0] not [y][0]...
        # v line -- actually worked before but only because LAMBDA2 === 1 so no epsilon dependence
        # v_pos = i*v_sgn*LAMBDA2(Lambda[0], Data[i-1][0])*STEP
        # v1_3d_pos = i*v_sgn*LAMBDA2(Lambda[1], Data[i-1][0])*STEP
        # v2_3d_pos = i*v_sgn*LAMBDA2(Lambda[2], Data[i-1][0])*STEP
        # Data[i][0] = DataPoint(Vec2D(0.0, v_pos), 0, 0, 1, 0,
        #                        [Vec2D(U0, v_pos + V0),
        #                         Vec2D(U0, v1_3d_pos + V0),
        #                         Vec2D(U0, v2_3d_pos + V0)],
        #                        [Vec2D(1.0, 0), Vec2D(1.0, 0), Vec2D(1.0, 0)],
        #                        [mathobj.KGaussian(U0, v_pos + V0, 0),
        #                         mathobj.KGaussian(U0, v1_3d_pos + V0, 1),
        #                         mathobj.KGaussian(U0, v2_3d_pos + V0, 2)],
        #                        0.0, 0.0, 0.0, 0.0)

        print(quadrant, v/PHISIZE)  # gives percent complete
        # at each u,v do a perp step, an actuated perp step, and a parallel step,
        # completing first an entire u row of perp/actuated steps, then looping
        # back and doing parallel steps, essentially advancing us one u line
        for u in range(0, THETASIZE):
            if hit_nan(Data[v][u]) or (u > 0 and hit_nan(Data[v][u-1])) or hit_nan(Data[v-1][u]):
                print("Data returned")
                # remove_v_singularities(Data, quadrant)
                return Data

            # check for singularity before proceeding:
            if at_perp_singularity(Data, u, v, v_sgn):
                print("hit singularity at v,u=", v, u)
                break

            # TODO NOTE new code here:
            # additional singularity check on epsilon/q?
            # if v > 2 and np.sign(Data[v-1][u].q) != np.sign(Data[v-2][u].q):
            #     print("dv epsilon changed signs...")
            #     return Data

            Data[v][u].epsilon = Data[v-1][u].epsilon + \
                v_sgn*STEP*(Data[v-1][u].q * Data[v-1][u].beta)

            # first compute d, then ensure d invertible
            try:
                perp_d = mathobj.d(Data[v-1][u], Data[v][u])
                # print(perp_d)
            except Exception as e:
                print(e)
                print("d not able to be calculated, ",
                      "likely M not invertible, i.e. hit singularity")
                # Data[v][u].epsilon = 1
                # remove_v_singularities(Data, quadrant)
                return Data

            # perp_d = mathobj.d(Data[v-1][u], Data[v][u])

            # now does both perp and actuated perp
            # TODO later set lambda somehwere globally visible...
            list_n_perp_step_double(
                mathobj, Data, u, v, v_sgn, perp_d, [0, 1, 2])

            if u != 0:
                if hit_nan(Data[v][u]) or (u > 0 and hit_nan(Data[v][u-1])) or hit_nan(Data[v-1][u]):
                    # Data[v][u].epsilon = 1
                    print("hit NaN at v,u=", v, u)
                    print("Data returned")
                    # remove_v_singularities(Data, quadrant)
                    return Data
                if at_parallel_singularity(Data, u, v, u_sgn):
                    print("hit singularity at v,u=", v, u)
                    # Data[v][u].epsilon = 1
                    break
                # parallel_d = mathobj.d(Data[v][u-1], Data[v][u]) # calculate separately so can start u at 0 and do dv along init v

                try:
                    parallel_d = mathobj.d(Data[v][u-1], Data[v][u])
                except Exception as e:
                    print(e)
                    print(
                        "d not able to be calculated, ",
                        "likely M not invertible, i.e. hit singularity")
                    # remove_v_singularities(Data, quadrant)
                    # Data[v][u].epsilon = 1
                    return Data

                n_step_double(Data, u, v, u_sgn, parallel_d)

    # finally remove v singularities just as with u above
    remove_v_singularities(Data, quadrant)

    # Done!
    print("Data returned")
    return Data


def remove_v_singularities(Data, quadrant):
    v_sgn = 1 if (quadrant == 1 or quadrant == 2) else -1
    v_singularity = 0
    for i in range(1, PHISIZE):
        if at_perp_singularity(Data, 0, i, v_sgn):
            print("hit singularity at v,u=", v, u)
            v_singularity = i
            break
    if v_singularity != 0:
        # if hit singularity, zero out past it:
        for i in range(v_singularity, PHISIZE):
            Data[i][0] = DataPoint(Vec2D(0, 0), 0.0, 0.0, 0.0, 0.0, [Vec2D(0.0, 0.0), Vec2D(0, 0), Vec2D(0, 0)],
                                   [Vec2D(0.0, 0.0), Vec2D(0, 0), Vec2D(0, 0)], [1.0, 1.0, 1.0], 0.0, 0.0, 0.0, 0.0)


def save_object(obj, filename, override):
    if override:
        with open(filename, 'wb') as output:
            # a - append
            # wb - override
            pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
    else:
        with open(filename, 'a') as output:
            # a - append
            # wb - override
            pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)


def save_state(Data):
    # for i in range(len(Data)):
    # try whole list of lists... if doesn't work try row by row
    save_object(Data, "Data.pkl", True)


def save_state2(Data):
    # for i in range(len(Data)):
    # try whole list of lists... if doesn't work try row by row
    # save_object(Data, "Data.pkl")
    for y in range(len(Data)):
        for x in range(len(Data[0])):
            save_object(Data[y][x], "Data.pkl", False)


def load_state(Data):
    with open("Data.pkl", "rb") as inpt:
        Data = pickle.load(inpt)
