# this file defines the relevant mathematics, both symbolically and then
# providing functions to extract the value at a point

from sympy import *
from consts import *
import consts
import math

u, v, alpha, beta, phi, theta, t, s = symbols(
    'u v alpha beta phi theta t s')

# additional symbols used just for calculating M^-1 and d:
# m components
m11, m12, m13, m21, m22, m23, m31, m32, m33 = symbols(
    'm11 m12 m13 m21 m22 m23 m31 m31 m33')
# k components - numeric value of curvature at pt
k1, k2, k3 = symbols('k1 k2 k3')
# variables used in Kbar
kbar1, kbar2, kbar3 = symbols('kbar1 kbar2 kbar3')


def construct_d_symb():
    # purely symbolic, at point plug in lambda_1,2, dup etc. and extract 0->1/betadvb 1->1/alphadus etc.
    M_mat = Matrix([[m11, m12, m13], [m21, m22, m23], [m31, m32, m33]])
    M_mat_inv = M_mat.inv()
    K_mat = Matrix([k1, k2, k3])
    Kbar_mat = Matrix([kbar1, kbar2, kbar3])
    # # try inv as:
    # try:
    #     M_mat.inv()
    # except:
    #     # handling code
    return M_mat_inv*(K_mat-Kbar_mat)


def construct_christoffel(g, i):
    # constructs the christoffel symbol T^i for the metric g (T as in Gamma ~ looks similar)
    # extract the terms of the first fundamental form from the metric:
    E = g[0, 0]
    G = g[1, 1]
    F = g[0, 1]
    # compute the desired symbol, calculating each tensor/matrix element then combining to output:
    if i == 1:  # i.e. for u/theta T^1
        c1_11 = (G*E.diff(theta) - 2*F*F.diff(theta) +
                 F*E.diff(phi))/(2*(E*G - F*F))
        c1_12 = (G*E.diff(phi) - F*G.diff(theta))/(2*(E*G-F*F))
        c1_22 = (2*G*F.diff(phi) - G*G.diff(theta) -
                 F*G.diff(phi)) / (2*(E*G-F*F))
        return (Matrix([[c1_11, c1_12], [c1_12, c1_22]]))
    else:  # i.e. T^2
        c2_11 = (2*E*F.diff(theta) - E*E.diff(phi) -
                 F*E.diff(theta))/(2*(E*G-F*F))
        c2_12 = (E*G.diff(theta) - F*E.diff(phi))/(2*(E*G-F*F))
        c2_22 = (E*G.diff(phi) - 2*F*F.diff(phi) +
                 F*G.diff(theta))/(2*(E*G-F*F))
        return (Matrix([[c2_11, c2_12], [c2_12, c2_22]]))


def construct_geodesic_curvature(mathobj, ut, vt, g, christoffel1, christoffel2):
    # returns an equation for calculating the geodesic curvature
    # of the line (u(t),v(t)) on surface of metric g
    E = g[0, 0]
    G = g[1, 1]
    F = g[0, 1]
    # see wolfram mathworld "Geodesic Curvature" page for following equation:
    out = (sqrt(E*G - F**2)*(-christoffel2[0, 0]*(Pow(ut.diff(t), 3)) +
                             christoffel1[1, 1]*(vt.diff(t)**3) -
                             (2*christoffel2[0, 1] - christoffel1[0, 0])*(Pow(ut.diff(t), 2))*vt.diff(t) +
                             (2*christoffel1[0, 1] - christoffel2[1, 1])*ut.diff(t)*(vt.diff(t)**2) +
                             ut.diff(t).diff(t)*vt.diff(t) -
                             vt.diff(t).diff(t)*ut.diff(t)) *
           ((E*ut.diff(t)**2 + 2*F*ut.diff(t)*vt.diff(t) +
             G*vt.diff(t)**2)**(-3/2)))
    return out


def construct_KGaussian(mathobj, surface):
    g = mathobj.g[surface]
    G2 = mathobj.christoffel2_symb[surface]
    gdescr = g[0, 0]*g[1, 1] - g[0, 1]*g[0, 1]
    return (1/sqrt(gdescr))*(((sqrt(gdescr)/g[0, 0])*G2[0, 0]).diff(phi) -
                             ((sqrt(gdescr)/g[0, 0])*G2[0, 1]).diff(theta))
    # ((sqrt(gdescr)/g[0, 0])*G2[0, 1]).diff(u)) NOTE TODO !!! this could have been root of errors


class MathObj:
    # this class merely bundles together all the needed formulae
    # and provides functions to get their values at a point
    def __init__(self, RS, curve1u, curve1v, curve2u, curve2v):
        # the mapping for the starting surface:
        self.Rs = RS  # plural R... else conflicts with function of same name :/
        # calculating the metric for all surfaces:
        self.g = []
        for i in range(len(RS)):
            self.g.append(Matrix([[(RS[i].diff(theta)).dot(RS[i].diff(theta)), (RS[i].diff(phi)).dot(RS[i].diff(theta))],
                                  [(RS[i].diff(phi)).dot(RS[i].diff(theta)), (RS[i].diff(phi)).dot(RS[i].diff(phi))]]))

        # calculating the terms of the first fundamental form for all surfaces:
        self.IE_symb = []
        for i in range(len(self.g)):
            self.IE_symb.append(RS[i].diff(theta).dot(RS[i].diff(theta)))

        # computing christoffel symbols:
        # split into its two subscriptable elements
        self.christoffel1_symb = []
        for i in range(len(self.g)):
            self.christoffel1_symb.append(construct_christoffel(self.g[i], 1))
        self.christoffel2_symb = []
        for i in range(len(self.g)):
            self.christoffel2_symb.append(construct_christoffel(self.g[i], 2))

        # generating equations for geodesic curvature along inital u and v curves:
        # TODO right?
        self.geodesic_curvature_u_symb = construct_geodesic_curvature(
            self, curve1u, curve1v, self.g[1], self.christoffel1_symb[1], self.christoffel2_symb[1])
        self.geodesic_curvature_v_symb = construct_geodesic_curvature(
            self, curve2u, curve2v, self.g[1], self.christoffel1_symb[1], self.christoffel2_symb[1])
        self.curve1u = curve1u
        self.curve1v = curve1v
        self.curve2u = curve2u
        self.curve2v = curve2v

        # now calculate gaussian curvature completely intrinsically
        # (i.e. using just the first fundamental form):
        self.KGaussian_symb = []
        for i in range(len(self.g)):
            self.KGaussian_symb.append(construct_KGaussian(self, i))

        # NOTE new
        self.d_symb = construct_d_symb()
        # --- end constructor ---

    def R(self, x, y, surface):
        return self.Rs[surface].subs([(theta, x), (phi, y)])

    # update this, should be a func of t not theta phi

    def geodesic_curvature(self, oldpoint, i, dist):
        # TODO fix for multiple surfaces
        if i == 1:
            return self.geodesic_curvature_u_symb.subs([(t, dist), (theta, oldpoint.r[1].x), (phi, oldpoint.r[1].y)])
        else:
            return self.geodesic_curvature_v_symb.subs([(t, dist), (theta, oldpoint.r[1].x), (phi, oldpoint.r[1].y)])

    def christoffel(self, oldpoint, surface, i):
        # for surface and the ith christoffel (i.e. 1 or 2)
        if i == 1:  # i.e. for u
            return self.christoffel1_symb[surface].subs([(theta, oldpoint.r[surface].x), (phi, oldpoint.r[surface].y)])
        else:
            return self.christoffel2_symb[surface].subs([(theta, oldpoint.r[surface].x), (phi, oldpoint.r[surface].y)])

    def KGaussian(self, uin, vin, surface):
        return self.KGaussian_symb[surface].subs([(theta, uin), (phi, vin)])

    def calc_log_deriv_term(self, oldpoint, new_point, big_lambda, small_lambda):
        # TODO NOTE is this correct for computing the necessary log lambda derivs?
        # TODO CHANGE set power in denominator to 1 not 2?
        if small_lambda == 1:
            # print(LAMBDA1(big_lambda, new_point),
            #       new_point.epsilon, oldpoint.q)
            d_top = math.log(LAMBDA1(big_lambda, new_point)) - \
                math.log(LAMBDA1(big_lambda, oldpoint))
            d_bottom = (LAMBDA2(big_lambda, oldpoint)**1) * \
                (new_point.epsilon - oldpoint.epsilon)
            # print(d_top)
            return -1*d_top/d_bottom
        else:
            d_top = math.log(LAMBDA2(big_lambda, new_point)) - \
                math.log(LAMBDA2(big_lambda, oldpoint))
            d_bottom = (LAMBDA1(big_lambda, oldpoint)**1) * \
                (new_point.epsilon - oldpoint.epsilon)
            # print(d_top)
            return -1*d_top/d_bottom

    def first_deriv_term(self, olderpoint, oldpoint, new_point, big_lambda):
        # calculate terms at new and old epsilon - since second derivative need an even older point if possible...
        # deriv1 = (LAMBDA2(big_lambda, new_point) - 2*LAMBDA2(big_lambda, oldpoint) + LAMBDA2(big_lambda,
        #                                                                                      olderpoint))/((new_point.epsilon-oldpoint.epsilon)*(oldpoint.epsilon - olderpoint.epsilon))
        # term1 = (1/LAMBDA1(big_lambda, oldpoint))*
        # will take shortcut since lambda2 not a func of epsilon -> dlambda2/depsilon = 0
        return 0

    def second_deriv_term(self, olderpoint, oldpoint, new_point, big_lambda):
        # again lambda2 not epsilon func, so becomes just 1/lambda2*secondderivlambda2
        deriv1 = (LAMBDA1(big_lambda, new_point) - 2*LAMBDA1(big_lambda, oldpoint) + LAMBDA1(big_lambda,
                                                                                             olderpoint))/((new_point.epsilon - oldpoint.epsilon)*(oldpoint.epsilon - olderpoint.epsilon))
        return deriv1  # since 1/lambda2 = 1/1 = 1

    def third_deriv_term(self, olderpoint, oldpoint, new_point, big_lambda):
        top = log((LAMBDA1(big_lambda, new_point) **
                   2)/LAMBDA2(big_lambda, new_point)) - log((LAMBDA1(big_lambda, oldpoint) **
                                                             2)/LAMBDA2(big_lambda, oldpoint))
        return top/((LAMBDA2(big_lambda, oldpoint)**2)*(new_point.epsilon - oldpoint.epsilon))

    def fourth_deriv_term(self, olderpoint, oldpoint, new_point, big_lambda):
        top = log((LAMBDA1(big_lambda, new_point)
                   )/(LAMBDA2(big_lambda, new_point)**2)) - log((LAMBDA1(big_lambda, oldpoint))/(LAMBDA2(big_lambda, oldpoint)**2))
        return top/((LAMBDA1(big_lambda, oldpoint)**2)*(new_point.epsilon - oldpoint.epsilon))

    def new_generate_kbar_term(self, olderpoint, oldpoint, new_point, big_lambda):
        term = -(oldpoint.b/LAMBDA2(big_lambda, oldpoint))**2
        term -= (oldpoint.s/LAMBDA1(big_lambda, oldpoint))**2
        term -= (1/(LAMBDA1(big_lambda, oldpoint) *
                    LAMBDA2(big_lambda, oldpoint)))*(oldpoint.p**2)*self.first_deriv_term(olderpoint, oldpoint, new_point, big_lambda)  # *DERIVSTUFF
        term -= (1/(LAMBDA1(big_lambda, oldpoint) *
                    LAMBDA2(big_lambda, oldpoint)))*(oldpoint.q**2)*self.second_deriv_term(olderpoint, oldpoint, new_point, big_lambda)  # *DERIVSTUFF
        term += oldpoint.b*oldpoint.q * \
            self.third_deriv_term(olderpoint, oldpoint, new_point,
                                  big_lambda)  # *DERIVSTUFF
        term += oldpoint.s*oldpoint.p * \
            self.fourth_deriv_term(olderpoint, oldpoint, new_point,
                                   big_lambda)  # *DERIVSTUFF
        return term

    def generate_kbar_term(self, oldpoint, new_point, big_lambda):
        # need du p is this an error comparing to page 3 u data?
        # so... TODO NOTE do we need du_p or dv_p... for now assuming du
        # fix term1, should take derivative e.g.... already did with dup but add lambda1 derivative?
        # term1 = -1*self.calc_log_deriv_term(oldpoint, new_point,
        #                                     big_lambda, 2)*(1/oldpoint.alpha)*(oldpoint.dup)
        # extra for calculating lambda2 deriv involving q:
        term0 = (oldpoint.b/(oldpoint.beta *
                             (LAMBDA2(big_lambda, oldpoint)**2)))
        term0 *= (oldpoint.q*(LAMBDA2(big_lambda, new_point)-LAMBDA2(big_lambda,
                                                                     oldpoint))/(new_point.epsilon - oldpoint.epsilon))
        term0 *= -1

        term1 = -1*self.calc_log_deriv_term(oldpoint, new_point,
                                            big_lambda, 2)*(1/oldpoint.alpha)*(oldpoint.dup)
        # do u deriv of lambda rather than epsilon, simpler
        # print("delta u: ",  new_point.uv.x - oldpoint.uv.x,
        #       " x y", new_point.uv.x, " ", new_point.uv.y, "oldpoint :", oldpoint.uv.x, " ", oldpoint.uv.y)
        # term1 += (oldpoint.s/(oldpoint.alpha *
        #                       (LAMBDA1(big_lambda, oldpoint)**2))) * ((LAMBDA1(big_lambda, new_point)-LAMBDA1(big_lambda, oldpoint))/(new_point.uv.x - oldpoint.uv.x))  # need du step size just = newpoint.u - oldpoint.u
        # doing d epsilon derivative:
        termadd = (oldpoint.s/(oldpoint.alpha *
                               (LAMBDA1(big_lambda, oldpoint)**2)))
        termaddmult = (oldpoint.p*(LAMBDA1(big_lambda, new_point)-LAMBDA1(big_lambda,
                                                                          oldpoint))/(new_point.epsilon - oldpoint.epsilon))
        # print(termadd*termaddmult)
        term1 += termadd*termaddmult
        term2 = -1*(oldpoint.b/LAMBDA2(big_lambda, oldpoint) +
                    self.calc_log_deriv_term(oldpoint, new_point, big_lambda, 1)*oldpoint.q)**2
        term3 = -1*(oldpoint.s/LAMBDA1(big_lambda, oldpoint) +
                    self.calc_log_deriv_term(oldpoint, new_point, big_lambda, 2)*oldpoint.p)**2
        # print(term1)
        return term0 + term1 + term2 + term3  # ensure all terms correct

    # NOTE new have to subs in a *lot* to make the equation work
    def d(self, oldpoint, new_point):
        # compute M terms:
        # 1/(LAMBDA2(consts.Lambda[0], oldpoint)**2)
        m11_val = LAMBDA2(consts.Lambda[0], oldpoint)**2
        # -1/(LAMBDA1(consts.Lambda[0], oldpoint)**2)
        m12_val = LAMBDA1(consts.Lambda[0], oldpoint)**2
        # TODO trickier, ensure works properly
        m13_val = self.calc_log_deriv_term(
            oldpoint, new_point, consts.Lambda[0], 1)
        # 1/(LAMBDA2(consts.Lambda[1], oldpoint)**2)
        m21_val = LAMBDA2(consts.Lambda[1], oldpoint)**2
        # -1/(LAMBDA1(consts.Lambda[1], oldpoint)**2)
        m22_val = LAMBDA1(consts.Lambda[1], oldpoint)**2
        m23_val = self.calc_log_deriv_term(
            oldpoint, new_point, consts.Lambda[1], 1)
        # 1/(LAMBDA2(consts.Lambda[1], oldpoint)**2)
        m31_val = LAMBDA2(consts.Lambda[2], oldpoint)**2
        # -1/(LAMBDA1(consts.Lambda[2], oldpoint)**2)
        m32_val = LAMBDA1(consts.Lambda[2], oldpoint)**2
        m33_val = self.calc_log_deriv_term(
            oldpoint, new_point, consts.Lambda[2], 1)
        # print(m11_val, " ", m12_val, " ", m13_val, " ", m21_val, " ", m22_val,
        #       " ", m23_val, " ", m31_val, " ", m32_val, " ", m33_val, " ")  # not the problem
        # compute K terms:
        k1_val = oldpoint.K[0]  # make sure terms right
        k2_val = oldpoint.K[1]  # TODO can generalize to higher dimensions...
        k3_val = oldpoint.K[2]
        # print(k3_val) # not the problem
        # compute Kbar terms:
        kbar1_val = self.generate_kbar_term(
            oldpoint, new_point, consts.Lambda[0])
        kbar2_val = self.generate_kbar_term(
            oldpoint, new_point, consts.Lambda[1])
        kbar3_val = self.generate_kbar_term(
            oldpoint, new_point, consts.Lambda[2])
        # print(kbar1_val, " ", kbar2_val, " ", kbar3_val)  # not the problem
        # d = M^-1.(K-Kbar)
        d_numeric = self.d_symb.subs([(m11, m11_val), (m12, m12_val), (m13, m13_val),
                                      (m21, m21_val), (m22,
                                                       m22_val), (m23, m23_val),
                                      (m31, m31_val), (m32,
                                                       m32_val), (m33, m33_val),
                                      (k1, k1_val), (k2, k2_val), (k3, k3_val),
                                      (kbar1, kbar1_val), (kbar2, kbar2_val), (kbar3, kbar3_val)])
        # print(d_numeric)  # the third entry is nan...
        return d_numeric

    def new_d(self, olderpoint, oldpoint, new_point):
        # compute M terms:
        # 1/(LAMBDA2(consts.Lambda[0], oldpoint)**2)
        m11_val = 1/(LAMBDA2(consts.Lambda[0], oldpoint)**2)
        # -1/(LAMBDA1(consts.Lambda[0], oldpoint)**2)
        m12_val = -1/(LAMBDA1(consts.Lambda[0], oldpoint)**2)
        # TODO trickier, ensure works properly
        m13_val = (-1/(LAMBDA1(consts.Lambda[0], oldpoint) * LAMBDA2(consts.Lambda[0], oldpoint)))*(
            (LAMBDA1(consts.Lambda[0], new_point)-LAMBDA1(consts.Lambda[0], oldpoint))/(LAMBDA2(consts.Lambda[0], oldpoint)*(new_point.epsilon-oldpoint.epsilon)))
        # self.calc_log_deriv_term(
        #     oldpoint, new_point, consts.Lambda[0], 1)
        # 1/(LAMBDA2(consts.Lambda[1], oldpoint)**2)
        m21_val = 1/LAMBDA2(consts.Lambda[1], oldpoint)**2
        # -1/(LAMBDA1(consts.Lambda[1], oldpoint)**2)
        m22_val = -1/LAMBDA1(consts.Lambda[1], oldpoint)**2
        m23_val = (-1/(LAMBDA1(consts.Lambda[1], oldpoint) * LAMBDA2(consts.Lambda[1], oldpoint)))*(
            (LAMBDA1(consts.Lambda[1], new_point)-LAMBDA1(consts.Lambda[1], oldpoint))/(LAMBDA2(consts.Lambda[1], oldpoint)*(new_point.epsilon-oldpoint.epsilon)))
        # 1/(LAMBDA2(consts.Lambda[1], oldpoint)**2)
        m31_val = 1/LAMBDA2(consts.Lambda[2], oldpoint)**2
        # -1/(LAMBDA1(consts.Lambda[2], oldpoint)**2)
        m32_val = -1/LAMBDA1(consts.Lambda[2], oldpoint)**2
        m33_val = (-1/(LAMBDA1(consts.Lambda[2], oldpoint) * LAMBDA2(consts.Lambda[2], oldpoint)))*(
            (LAMBDA1(consts.Lambda[2], new_point)-LAMBDA1(consts.Lambda[2], oldpoint))/(LAMBDA2(consts.Lambda[2], oldpoint)*(new_point.epsilon-oldpoint.epsilon)))
        # print(m11_val, " ", m12_val, " ", m13_val, " ", m21_val, " ", m22_val,
        #       " ", m23_val, " ", m31_val, " ", m32_val, " ", m33_val, " ")  # not the problem
        # compute K terms:
        k1_val = oldpoint.K[0]  # make sure terms right
        k2_val = oldpoint.K[1]  # TODO can generalize to higher dimensions...
        k3_val = oldpoint.K[2]
        # print(k3_val) # not the problem
        # compute Kbar terms:
        kbar1_val = self.new_generate_kbar_term(
            olderpoint, oldpoint, new_point, consts.Lambda[0])
        kbar2_val = self.new_generate_kbar_term(
            olderpoint, oldpoint, new_point, consts.Lambda[1])
        kbar3_val = self.new_generate_kbar_term(
            olderpoint, oldpoint, new_point, consts.Lambda[2])
        # print(kbar1_val, " ", kbar2_val, " ", kbar3_val)  # not the problem
        # d = M^-1.(K-Kbar)
        d_numeric = self.d_symb.subs([(m11, m11_val), (m12, m12_val), (m13, m13_val),
                                      (m21, m21_val), (m22,
                                                       m22_val), (m23, m23_val),
                                      (m31, m31_val), (m32,
                                                       m32_val), (m33, m33_val),
                                      (k1, k1_val), (k2, k2_val), (k3, k3_val),
                                      (kbar1, kbar1_val), (kbar2, kbar2_val), (kbar3, kbar3_val)])
        # print(d_numeric)
        return d_numeric
