from consts import *
from mayavi import mlab
import json

# removed old plotting functions, now plotting all surfaces along with epsilon coloring


def plot_3d_4_quadrants_2_surfaces_epsilon(mathobj, Data1, Data2, Data3, Data4):
    # repeat the following for all 4 data sets:
    epsilon1_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            epsilon1_data[y][x] = float(Data1[y][x].epsilon)
            if epsilon1_data[y][x] == 1:
                epsilon1_data[y][x] = np.nan
    epsilon2_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            epsilon2_data[y][x] = float(Data2[y][x].epsilon)
            if epsilon2_data[y][x] == 1:
                epsilon2_data[y][x] = np.nan
    epsilon3_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            epsilon3_data[y][x] = float(Data3[y][x].epsilon)
            if epsilon3_data[y][x] == 1:
                epsilon3_data[y][x] = np.nan
    epsilon4_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            epsilon4_data[y][x] = float(Data4[y][x].epsilon)
            if epsilon4_data[y][x] == 1:
                epsilon4_data[y][x] = np.nan

    epsilon1_data = np.asarray(epsilon1_data)
    epsilon2_data = np.asarray(epsilon2_data)
    epsilon3_data = np.asarray(epsilon3_data)
    epsilon4_data = np.asarray(epsilon4_data)

    R1_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R1_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R1_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            R1_x_data[y][x] = float(mathobj.R(Data1[y][x].r[1].x, Data1[y][x].r[1].y, 1)[
                0])
            R1_y_data[y][x] = float(mathobj.R(Data1[y][x].r[1].x,
                                              Data1[y][x].r[1].y, 1)[1])
            R1_z_data[y][x] = float(mathobj.R(Data1[y][x].r[1].x,
                                              Data1[y][x].r[1].y, 1)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    R1_x_out = np.asarray(R1_x_data)
    R1_y_out = np.asarray(R1_y_data)
    R1_z_out = np.asarray(R1_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if R1_x_out[y][x] == 0 and R1_y_out[y][x] == 0 and R1_z_out[y][x] == 1 and (x != 0 or y != 0):
                R1_x_out[y][x] = np.nan
                R1_y_out[y][x] = np.nan
                R1_z_out[y][x] = np.nan
    # again for second surface
    R1_x_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R1_y_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R1_z_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            R1_x_data2[y][x] = float(mathobj.R(Data1[y][x].r[2].x, Data1[y][x].r[2].y, 2)[
                0])
            R1_y_data2[y][x] = float(mathobj.R(Data1[y][x].r[2].x,
                                               Data1[y][x].r[2].y, 2)[1])
            R1_z_data2[y][x] = float(mathobj.R(Data1[y][x].r[2].x,
                                               Data1[y][x].r[2].y, 2)[2]) + 0.5

            # R1_x_data2[y][x] = float(mathobj.RM2(Data1[y][x].r2omega.x, Data1[y][x].r2omega.y)[
            #     0])
            # R1_y_data2[y][x] = float(mathobj.RM2(Data1[y][x].r2omega.x,
            #                              Data1[y][x].r2omega.y)[1])
            # R1_z_data2[y][x] = float(mathobj.RM2(Data1[y][x].r2omega.x,
            #                              Data1[y][x].r2omega.y)[2])
            # print(R1_x_data[y][x], R1_y_data[y][x], R1_z_data[y][x])
    # convert to numpy array - just since matplotlib expects this format to plot
    R1_x_out2 = np.asarray(R1_x_data2)
    R1_y_out2 = np.asarray(R1_y_data2)
    R1_z_out2 = np.asarray(R1_z_data2)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if R1_x_out2[y][x] == 1 and R1_y_out2[y][x] == 0 and R1_z_out2[y][x] == 0 and (x != 0 or y != 0):
                R1_x_out2[y][x] = np.nan
                R1_y_out2[y][x] = np.nan
                R1_z_out2[y][x] = np.nan

    # now extract 3d initial curves:
    RS1_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS1_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS1_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            RS1_x_data[y][x] = float(mathobj.R(Data1[y][x].r[0].x, Data1[y][x].r[0].y, 0)[
                0])
            RS1_y_data[y][x] = float(mathobj.R(Data1[y][x].r[0].x,
                                               Data1[y][x].r[0].y, 0)[1])
            RS1_z_data[y][x] = float(mathobj.R(Data1[y][x].r[0].x,
                                               Data1[y][x].r[0].y, 0)[2])

            # RS1_x_data[y][x] = float(mathobj.RS(Data1[y][x].r.x, Data1[y][x].r.y)[
            #     0])
            # RS1_y_data[y][x] = float(mathobj.RS(Data1[y][x].r.x,
            #                             Data1[y][x].r.y)[1])
            # RS1_z_data[y][x] = float(mathobj.RS(Data1[y][x].r.x,
            #                             Data1[y][x].r.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    RS1_x_out = np.asarray(RS1_x_data)
    RS1_y_out = np.asarray(RS1_y_data)
    RS1_z_out = np.asarray(RS1_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if RS1_x_out[y][x] == 0 and RS1_y_out[y][x] == 0 and RS1_z_out[y][x] == 0 and (x != 0 or y != 0):
                RS1_x_out[y][x] = np.nan
                RS1_y_out[y][x] = np.nan
                RS1_z_out[y][x] = np.nan

    # ------

    R2_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R2_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R2_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            R2_x_data[y][x] = float(mathobj.R(Data2[y][x].r[1].x, Data2[y][x].r[1].y, 1)[
                0])
            R2_y_data[y][x] = float(mathobj.R(Data2[y][x].r[1].x,
                                              Data2[y][x].r[1].y, 1)[1])
            R2_z_data[y][x] = float(mathobj.R(Data2[y][x].r[1].x,
                                              Data2[y][x].r[1].y, 1)[2])
            # R2_x_data[y][x] = float(mathobj.RM(Data2[y][x].romega.x, Data2[y][x].romega.y)[
            #     0])
            # R2_y_data[y][x] = float(mathobj.RM(Data2[y][x].romega.x,
            #                            Data2[y][x].romega.y)[1])
            # R2_z_data[y][x] = float(mathobj.RM(Data2[y][x].romega.x,
            #                            Data2[y][x].romega.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    R2_x_out = np.asarray(R2_x_data)
    R2_y_out = np.asarray(R2_y_data)
    R2_z_out = np.asarray(R2_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if R2_x_out[y][x] == 0 and R2_y_out[y][x] == 0 and R2_z_out[y][x] == 1 and (x != 0 or y != 0):
                R2_x_out[y][x] = np.nan
                R2_y_out[y][x] = np.nan
                R2_z_out[y][x] = np.nan

    # again for second surface
    R2_x_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R2_y_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R2_z_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            R2_x_data2[y][x] = float(mathobj.R(Data2[y][x].r[2].x, Data2[y][x].r[2].y, 2)[
                0])
            R2_y_data2[y][x] = float(mathobj.R(Data2[y][x].r[2].x,
                                               Data2[y][x].r[2].y, 2)[1])
            R2_z_data2[y][x] = float(mathobj.R(Data2[y][x].r[2].x,
                                               Data2[y][x].r[2].y, 2)[2]) + 0.5
            # R2_x_data2[y][x] = float(mathobj.RM2(Data2[y][x].r2omega.x, Data2[y][x].r2omega.y)[
            #     0])
            # R2_y_data2[y][x] = float(mathobj.RM2(Data2[y][x].r2omega.x,
            #                              Data2[y][x].r2omega.y)[1])
            # R2_z_data2[y][x] = float(mathobj.RM2(Data2[y][x].r2omega.x,
            #                              Data2[y][x].r2omega.y)[2])
            # print(R1_x_data[y][x], R1_y_data[y][x], R1_z_data[y][x])
    # convert to numpy array - just since matplotlib expects this format to plot
    R2_x_out2 = np.asarray(R2_x_data2)
    R2_y_out2 = np.asarray(R2_y_data2)
    R2_z_out2 = np.asarray(R2_z_data2)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if R2_x_out2[y][x] == 1 and R2_y_out2[y][x] == 0 and R2_z_out2[y][x] == 0 and (x != 0 or y != 0):
                R2_x_out2[y][x] = np.nan
                R2_y_out2[y][x] = np.nan
                R2_z_out2[y][x] = np.nan

    # now extract 3d initial curves:
    RS2_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS2_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS2_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            RS2_x_data[y][x] = float(mathobj.R(Data2[y][x].r[0].x, Data2[y][x].r[0].y, 0)[
                0])
            RS2_y_data[y][x] = float(mathobj.R(Data2[y][x].r[0].x,
                                               Data2[y][x].r[0].y, 0)[1])
            RS2_z_data[y][x] = float(mathobj.R(Data2[y][x].r[0].x,
                                               Data2[y][x].r[0].y, 0)[2])
            # RS2_x_data[y][x] = float(mathobj.RS(Data2[y][x].r.x, Data2[y][x].r.y)[
            #     0])
            # RS2_y_data[y][x] = float(mathobj.RS(Data2[y][x].r.x,
            #                             Data2[y][x].r.y)[1])
            # RS2_z_data[y][x] = float(mathobj.RS(Data2[y][x].r.x,
            #                             Data2[y][x].r.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    RS2_x_out = np.asarray(RS2_x_data)
    RS2_y_out = np.asarray(RS2_y_data)
    RS2_z_out = np.asarray(RS2_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if RS2_x_out[y][x] == 0 and RS2_y_out[y][x] == 0 and RS2_z_out[y][x] == 0 and (x != 0 or y != 0):
                RS2_x_out[y][x] = np.nan
                RS2_y_out[y][x] = np.nan
                RS2_z_out[y][x] = np.nan

    # ------

    R3_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R3_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R3_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            R3_x_data[y][x] = float(mathobj.R(Data3[y][x].r[1].x, Data3[y][x].r[1].y, 1)[
                0])
            R3_y_data[y][x] = float(mathobj.R(Data3[y][x].r[1].x,
                                              Data3[y][x].r[1].y, 1)[1])
            R3_z_data[y][x] = float(mathobj.R(Data3[y][x].r[1].x,
                                              Data3[y][x].r[1].y, 1)[2])

            # R3_x_data[y][x] = float(mathobj.RM(Data3[y][x].romega.x, Data3[y][x].romega.y)[
            #     0])
            # R3_y_data[y][x] = float(mathobj.RM(Data3[y][x].romega.x,
            #                            Data3[y][x].romega.y)[1])
            # R3_z_data[y][x] = float(mathobj.RM(Data3[y][x].romega.x,
            #                            Data3[y][x].romega.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    R3_x_out = np.asarray(R3_x_data)
    R3_y_out = np.asarray(R3_y_data)
    R3_z_out = np.asarray(R3_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if R3_x_out[y][x] == 0 and R3_y_out[y][x] == 0 and R3_z_out[y][x] == 1 and (x != 0 or y != 0):
                R3_x_out[y][x] = np.nan
                R3_y_out[y][x] = np.nan
                R3_z_out[y][x] = np.nan

    # again for second surface
    R3_x_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R3_y_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R3_z_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            R3_x_data2[y][x] = float(mathobj.R(Data3[y][x].r[2].x, Data3[y][x].r[2].y, 2)[
                0])
            R3_y_data2[y][x] = float(mathobj.R(Data3[y][x].r[2].x,
                                               Data3[y][x].r[2].y, 2)[1])
            R3_z_data2[y][x] = float(mathobj.R(Data3[y][x].r[2].x,
                                               Data3[y][x].r[2].y, 2)[2]) + 0.5

            # R3_x_data2[y][x] = float(mathobj.RM2(Data3[y][x].r2omega.x, Data3[y][x].r2omega.y)[
            #     0])
            # R3_y_data2[y][x] = float(mathobj.RM2(Data3[y][x].r2omega.x,
            #                              Data3[y][x].r2omega.y)[1])
            # R3_z_data2[y][x] = float(mathobj.RM2(Data3[y][x].r2omega.x,
            #                              Data3[y][x].r2omega.y)[2])
            # print(R1_x_data[y][x], R1_y_data[y][x], R1_z_data[y][x])
    # convert to numpy array - just since matplotlib expects this format to plot
    R3_x_out2 = np.asarray(R3_x_data2)
    R3_y_out2 = np.asarray(R3_y_data2)
    R3_z_out2 = np.asarray(R3_z_data2)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if R3_x_out2[y][x] == 1 and R3_y_out2[y][x] == 0 and R3_z_out2[y][x] == 0 and (x != 0 or y != 0):
                R3_x_out2[y][x] = np.nan
                R3_y_out2[y][x] = np.nan
                R3_z_out2[y][x] = np.nan

    # now extract 3d initial curves:
    RS3_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS3_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS3_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            RS3_x_data[y][x] = float(mathobj.R(Data3[y][x].r[0].x, Data3[y][x].r[0].y, 0)[
                0])
            RS3_y_data[y][x] = float(mathobj.R(Data3[y][x].r[0].x,
                                               Data3[y][x].r[0].y, 0)[1])
            RS3_z_data[y][x] = float(mathobj.R(Data3[y][x].r[0].x,
                                               Data3[y][x].r[0].y, 0)[2])

            # RS3_x_data[y][x] = float(mathobj.RS(Data3[y][x].r.x, Data3[y][x].r.y)[
            #     0])
            # RS3_y_data[y][x] = float(mathobj.RS(Data3[y][x].r.x,
            #                             Data3[y][x].r.y)[1])
            # RS3_z_data[y][x] = float(mathobj.RS(Data3[y][x].r.x,
            #                             Data3[y][x].r.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    RS3_x_out = np.asarray(RS3_x_data)
    RS3_y_out = np.asarray(RS3_y_data)
    RS3_z_out = np.asarray(RS3_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if RS3_x_out[y][x] == 0 and RS3_y_out[y][x] == 0 and RS3_z_out[y][x] == 0 and (x != 0 or y != 0):
                RS3_x_out[y][x] = np.nan
                RS3_y_out[y][x] = np.nan
                RS3_z_out[y][x] = np.nan

    # ------

    R4_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R4_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R4_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            R4_x_data[y][x] = float(mathobj.R(Data4[y][x].r[1].x, Data4[y][x].r[1].y, 1)[
                0])
            R4_y_data[y][x] = float(mathobj.R(Data4[y][x].r[1].x,
                                              Data4[y][x].r[1].y, 1)[1])
            R4_z_data[y][x] = float(mathobj.R(Data4[y][x].r[1].x,
                                              Data4[y][x].r[1].y, 1)[2])

            # R4_x_data[y][x] = float(mathobj.RM(Data4[y][x].romega.x, Data4[y][x].romega.y)[
            #     0])
            # R4_y_data[y][x] = float(mathobj.RM(Data4[y][x].romega.x,
            #                            Data4[y][x].romega.y)[1])
            # R4_z_data[y][x] = float(mathobj.RM(Data4[y][x].romega.x,
            #                            Data4[y][x].romega.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    R4_x_out = np.asarray(R4_x_data)
    R4_y_out = np.asarray(R4_y_data)
    R4_z_out = np.asarray(R4_z_data)
    # remove 0 elements from line plot(caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if R4_x_out[y][x] == 0 and R4_y_out[y][x] == 0 and R4_z_out[y][x] == 1 and (x != 0 or y != 0):
                R4_x_out[y][x] = np.nan
                R4_y_out[y][x] = np.nan
                R4_z_out[y][x] = np.nan

    # again for second surface
    R4_x_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R4_y_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R4_z_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            R4_x_data2[y][x] = float(mathobj.R(Data4[y][x].r[2].x, Data4[y][x].r[2].y, 2)[
                0])
            R4_y_data2[y][x] = float(mathobj.R(Data4[y][x].r[2].x,
                                               Data4[y][x].r[2].y, 2)[1])
            R4_z_data2[y][x] = float(mathobj.R(Data4[y][x].r[2].x,
                                               Data4[y][x].r[2].y, 2)[2]) + 0.5

            # R4_x_data2[y][x] = float(mathobj.RM2(Data4[y][x].r2omega.x, Data4[y][x].r2omega.y)[
            #     0])
            # R4_y_data2[y][x] = float(mathobj.RM2(Data4[y][x].r2omega.x,
            #                              Data4[y][x].r2omega.y)[1])
            # R4_z_data2[y][x] = float(mathobj.RM2(Data4[y][x].r2omega.x,
            #                              Data4[y][x].r2omega.y)[2])
            # print(R1_x_data[y][x], R1_y_data[y][x], R1_z_data[y][x])
    # convert to numpy array - just since matplotlib expects this format to plot
    R4_x_out2 = np.asarray(R4_x_data2)
    R4_y_out2 = np.asarray(R4_y_data2)
    R4_z_out2 = np.asarray(R4_z_data2)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if R4_x_out2[y][x] == 1 and R4_y_out2[y][x] == 0 and R4_z_out2[y][x] == 0 and (x != 0 or y != 0):
                R4_x_out2[y][x] = np.nan
                R4_y_out2[y][x] = np.nan
                R4_z_out2[y][x] = np.nan

    # now extract 3d initial curves:
    RS4_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS4_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS4_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            RS4_x_data[y][x] = float(mathobj.R(Data4[y][x].r[0].x, Data4[y][x].r[0].y, 0)[
                0])
            RS4_y_data[y][x] = float(mathobj.R(Data4[y][x].r[0].x,
                                               Data4[y][x].r[0].y, 0)[1])
            RS4_z_data[y][x] = float(mathobj.R(Data4[y][x].r[0].x,
                                               Data4[y][x].r[0].y, 0)[2])
            # RS4_x_data[y][x] = float(mathobj.RS(Data4[y][x].r.x, Data4[y][x].r.y)[
            #     0])
            # RS4_y_data[y][x] = float(mathobj.RS(Data4[y][x].r.x,
            #                             Data4[y][x].r.y)[1])
            # RS4_z_data[y][x] = float(mathobj.RS(Data4[y][x].r.x,
            #                             Data4[y][x].r.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    RS4_x_out = np.asarray(RS4_x_data)
    RS4_y_out = np.asarray(RS4_y_data)
    RS4_z_out = np.asarray(RS4_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if RS4_x_out[y][x] == 0 and RS4_y_out[y][x] == 0 and RS4_z_out[y][x] == 0 and (x != 0 or y != 0):
                RS4_x_out[y][x] = np.nan
                RS4_y_out[y][x] = np.nan
                RS4_z_out[y][x] = np.nan

    # --------

    vspace = np.linspace(-np.pi, np.pi, 100)
    uspace = np.linspace(-np.pi, np.pi, 100)
    xspace, yspace = np.meshgrid(uspace, vspace)
    zspace = xspace*yspace*0
    # what to do for sphere...? # fix bounds... see old version on box
    x2 = 1 * np.outer(np.cos(uspace), np.sin(vspace))
    y2 = 1 * np.outer(np.sin(vspace), np.sin(uspace))
    z2space = 1 + 1 * np.outer(np.ones(np.size(vspace)), np.cos(vspace))
    z3space = 2 + np.exp(-((xspace**2)/(4) + 2*yspace**2))
    # Plot the surface, both the initial and final
    mlab.figure(bgcolor=(1, 1, 1))
    # mlab.mesh(xspace, yspace, zspace, color=(
    #     1.0, 0.7, 1.0), opacity=0.8)  # initial
    # mlab.mesh(x2, y2, z2space, color=(
    #     0.8, 0.7, 0.4), opacity=0.8)  # final
    # now need a third surface too
    # mlab.mesh(xspace, yspace, z3space, color=(
    #     0.3, 0.4, 0.8), opacity=0.8)

    # plotting the u and v initial curves
    mlab.plot3d(R1_x_out[0], R1_y_out[0], R1_z_out[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(R1_x_out[:, 0], R1_y_out[:, 0],
                R1_z_out[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(R2_x_out[0], R2_y_out[0], R2_z_out[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(R2_x_out[:, 0], R2_y_out[:, 0],
                R2_z_out[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(R3_x_out[0], R2_y_out[0], R3_z_out[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(R3_x_out[:, 0], R2_y_out[:, 0],
                R3_z_out[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(R4_x_out[0], R4_y_out[0], R3_z_out[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(R4_x_out[:, 0], R4_y_out[:, 0],
                R4_z_out[:, 0], color=(1.0, 0.6, 0.2))

    mlab.plot3d(R1_x_out2[0], R1_y_out2[0],
                R1_z_out2[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(R1_x_out2[:, 0], R1_y_out2[:, 0],
                R1_z_out2[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(R2_x_out2[0], R2_y_out2[0],
                R2_z_out2[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(R2_x_out2[:, 0], R2_y_out2[:, 0],
                R2_z_out2[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(R3_x_out2[0], R2_y_out2[0],
                R3_z_out2[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(R3_x_out2[:, 0], R2_y_out2[:, 0],
                R3_z_out2[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(R4_x_out2[0], R4_y_out2[0],
                R3_z_out2[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(R4_x_out2[:, 0], R4_y_out2[:, 0],
                R4_z_out2[:, 0], color=(1.0, 0.6, 0.2))

    mlab.plot3d(RS1_x_out[0], RS1_y_out[0],
                RS1_z_out[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(RS1_x_out[:, 0], RS1_y_out[:, 0],
                RS1_z_out[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(RS2_x_out[0], RS2_y_out[0],
                RS2_z_out[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(RS2_x_out[:, 0], RS2_y_out[:, 0],
                RS2_z_out[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(RS3_x_out[0], RS3_y_out[0],
                RS3_z_out[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(RS3_x_out[:, 0], RS3_y_out[:, 0],
                RS3_z_out[:, 0], color=(1.0, 0.6, 0.2))
    mlab.plot3d(RS4_x_out[0], RS4_y_out[0],
                RS4_z_out[0], color=(1.0, 0.0, 0.0))
    mlab.plot3d(RS4_x_out[:, 0], RS4_y_out[:, 0],
                RS4_z_out[:, 0], color=(1.0, 0.6, 0.2))

    # try clamping?

    # # TODO annoying matplotlib wont accept [:, y] for phisize!=thetasize, so flip order of xy in arrays just here for printing...
    # for x in range(1, THETASIZE):
    #     for y in range(1, PHISIZE):
    #         R1_x_out[y][x] = R1_x_out[x][y]
    #         R1_y_out[y][x] = R1_y_out[x][y]
    #         R1_z_out[y][x] = R1_z_out[x][y]
    #         R1_x_out2[y][x] = R1_x_out2[x][y]
    #         R1_y_out2[y][x] = R1_y_out2[x][y]
    #         R1_z_out2[y][x] = R1_y_out2[x][y]
    #         RS1_x_out[y][x] = RS1_x_out[x][y]
    #         RS1_y_out[y][x] = RS1_y_out[x][y]
    #         RS1_z_out[y][x] = RS1_z_out[x][y]
    #         epsilon1_data[y][x] = epsilon1_data[x][y]
    #         R2_x_out[y][x] = R2_x_out[x][y]
    #         R2_y_out[y][x] = R2_y_out[x][y]
    #         R2_z_out[y][x] = R2_z_out[x][y]
    #         R2_x_out2[y][x] = R2_x_out2[x][y]
    #         R2_y_out2[y][x] = R2_y_out2[x][y]
    #         R2_z_out2[y][x] = R2_y_out2[x][y]
    #         RS2_x_out[y][x] = RS2_x_out[x][y]
    #         RS2_y_out[y][x] = RS2_y_out[x][y]
    #         RS2_z_out[y][x] = RS2_z_out[x][y]
    #         epsilon2_data[y][x] = epsilon2_data[x][y]
    #         R3_x_out[y][x] = R3_x_out[x][y]
    #         R3_y_out[y][x] = R3_y_out[x][y]
    #         R3_z_out[y][x] = R3_z_out[x][y]
    #         R3_x_out2[y][x] = R3_x_out2[x][y]
    #         R3_y_out2[y][x] = R3_y_out2[x][y]
    #         R3_z_out2[y][x] = R3_y_out2[x][y]
    #         RS3_x_out[y][x] = RS3_x_out[x][y]
    #         RS3_y_out[y][x] = RS3_y_out[x][y]
    #         RS3_z_out[y][x] = RS3_z_out[x][y]
    #         epsilon3_data[y][x] = epsilon3_data[x][y]
    #         R4_x_out[y][x] = R4_x_out[x][y]
    #         R4_y_out[y][x] = R4_y_out[x][y]
    #         R4_z_out[y][x] = R4_z_out[x][y]
    #         R4_x_out2[y][x] = R4_x_out2[x][y]
    #         R4_y_out2[y][x] = R4_y_out2[x][y]
    #         R4_z_out2[y][x] = R4_y_out2[x][y]
    #         RS4_x_out[y][x] = RS4_x_out[x][y]
    #         RS4_y_out[y][x] = RS4_y_out[x][y]
    #         RS4_z_out[y][x] = RS4_z_out[x][y]
    #         epsilon4_data[y][x] = epsilon4_data[x][y]

    # now plot actual data of director curves
    for x in range(1, THETASIZE, 1):
        # TODO NOTE flip so we print constant u, varying v to get proper behavior for baromorph
        # plotting constant u lines - was columns but wrong.. should be rows.. u lines implies *varying* u not constant u!
        # mlab.points3d(R1_x_out[y], R1_y_out[y], R1_z_out[y], epsilon1_data[y], colormap='autumn')
        # mlab.points3d(R2_x_out[y], R2_y_out[y], R2_z_out[y], epsilon2_data[y], colormap='autumn')
        # mlab.points3d(R3_x_out[y], R3_y_out[y], R3_z_out[y], epsilon3_data[y], colormap='autumn')
        # mlab.points3d(R4_x_out[y], R4_y_out[y], R4_z_out[y], epsilon4_data[y], colormap='autumn')

        # mlab.points3d(R1_x_out2[y], R1_y_out2[y], R1_z_out2[y], epsilon1_data[y], colormap='autumn')
        # mlab.points3d(R2_x_out2[y], R2_y_out2[y], R2_z_out2[y], epsilon2_data[y], colormap='autumn')
        # mlab.points3d(R3_x_out2[y], R3_y_out2[y], R3_z_out2[y], epsilon3_data[y], colormap='autumn')
        # mlab.points3d(R4_x_out2[y], R4_y_out2[y], R4_z_out2[y], epsilon4_data[y], colormap='autumn')

        # mlab.points3d(RS1_x_out[y], RS1_y_out[y], RS1_z_out[y], epsilon1_data[y], colormap='autumn')
        # mlab.points3d(RS2_x_out[y], RS2_y_out[y], RS2_z_out[y], epsilon2_data[y], colormap='autumn')
        # mlab.points3d(RS3_x_out[y], RS3_y_out[y], RS3_z_out[y], epsilon3_data[y], colormap='autumn')
        # mlab.points3d(RS4_x_out[y], RS4_y_out[y], RS4_z_out[y], epsilon4_data[y], colormap='autumn')
        mlab.plot3d(R1_x_out[:, x], R1_y_out[:, x], R1_z_out[:, x],
                    epsilon1_data[:, x], colormap='autumn')
        mlab.plot3d(R2_x_out[:, x], R2_y_out[:, x], R2_z_out[:, x],
                    epsilon2_data[:, x], colormap='autumn')
        mlab.plot3d(R3_x_out[:, x], R3_y_out[:, x], R3_z_out[:, x],
                    epsilon3_data[:, x], colormap='autumn')
        mlab.plot3d(R4_x_out[:, x], R4_y_out[:, x], R4_z_out[:, x],
                    epsilon4_data[:, x], colormap='autumn')

        mlab.plot3d(R1_x_out2[:, x], R1_y_out2[:, x], R1_z_out2[:, x],
                    epsilon1_data[:, x], colormap='autumn')
        mlab.plot3d(R2_x_out2[:, x], R2_y_out2[:, x], R2_z_out2[:, x],
                    epsilon2_data[:, x], colormap='autumn')
        mlab.plot3d(R3_x_out2[:, x], R3_y_out2[:, x], R3_z_out2[:, x],
                    epsilon3_data[:, x], colormap='autumn')
        mlab.plot3d(R4_x_out2[:, x], R4_y_out2[:, x], R4_z_out2[:, x],
                    epsilon4_data[:, x], colormap='autumn')

        mlab.plot3d(RS1_x_out[:, x], RS1_y_out[:, x], RS1_z_out[:, x],
                    epsilon1_data[:, x], colormap='autumn')
        mlab.plot3d(RS2_x_out[:, x], RS2_y_out[:, x], RS2_z_out[:, x],
                    epsilon2_data[:, x], colormap='autumn')
        mlab.plot3d(RS3_x_out[:, x], RS3_y_out[:, x], RS3_z_out[:, x],
                    epsilon3_data[:, x], colormap='autumn')
        mlab.plot3d(RS4_x_out[:, x], RS4_y_out[:, x], RS4_z_out[:, x],
                    epsilon4_data[:, x], colormap='autumn')

        # try plotting with colormap... extract epsilon at each point
    mlab.show()


def json_out(mathobj, Data):
    # NOTE updated to accomodate 3 surfaces total
    # now do not conver to numpy, json does not recognize...
    # first decompose to an array of r values as in plotting
    R_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            # we then store the x, y, and z coordinates of the final surface at each point:
            R_x_data[y][x] = float(mathobj.RM(Data[y][x].romega.x, Data[y][x].romega.y)[
                0])
            R_y_data[y][x] = float(mathobj.RM(Data[y][x].romega.x,
                                              Data[y][x].romega.y)[1])
            R_z_data[y][x] = float(mathobj.RM(Data[y][x].romega.x,
                                              Data[y][x].romega.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot:
    R_x_out = np.asarray(R_x_data)
    R_y_out = np.asarray(R_y_data)
    R_z_out = np.asarray(R_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            # SEE HERE: this must be tweaked for the particular surface and initial point
            # e.g. the paraboloid mapping goes to (0,0,0) at the initial point (OR data just zeroed out), so that is what
            # must be checked for if we are not at the initial point
            if R_x_out[y][x] == 0 and R_y_out[y][x] == 0 and R_z_out[y][x] == 0 and x != 0 and y != 0:
                # np.nan gets the plotter to ignore this data without messing up array size
                R_x_out[y][x] = np.nan
                R_y_out[y][x] = np.nan
                R_z_out[y][x] = np.nan

    R_x_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R_y_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R_z_data2 = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            # we then store the x, y, and z coordinates of the final surface at each point:
            R_x_data2[y][x] = float(mathobj.RM2(Data[y][x].r2omega.x, Data[y][x].r2omega.y)[
                0])
            R_y_data2[y][x] = float(mathobj.RM2(Data[y][x].r2omega.x,
                                                Data[y][x].r2omega.y)[1])
            R_z_data2[y][x] = float(mathobj.RM2(Data[y][x].r2omega.x,
                                                Data[y][x].r2omega.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot:
    R_x_out2 = np.asarray(R_x_data2)
    R_y_out2 = np.asarray(R_y_data2)
    R_z_out2 = np.asarray(R_z_data2)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            # SEE HERE: this must be tweaked for the particular surface and initial point
            # e.g. the paraboloid mapping goes to (0,0,0) at the initial point (OR data just zeroed out), so that is what
            # must be checked for if we are not at the initial point
            if R_x_out2[y][x] == 0 and R_y_out2[y][x] == 0 and R_z_out2[y][x] == 0 and x != 0 and y != 0:
                # np.nan gets the plotter to ignore this data without messing up array size
                R_x_out2[y][x] = np.nan
                R_y_out2[y][x] = np.nan
                R_z_out2[y][x] = np.nan

    # now extract 3d initial curves (vs above is for final curves):
    RS_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS_z_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            RS_x_data[y][x] = float(mathobj.RS(Data[y][x].r.x, Data[y][x].r.y)[
                0])
            RS_y_data[y][x] = float(mathobj.RS(Data[y][x].r.x,
                                               Data[y][x].r.y)[1])
            RS_z_data[y][x] = float(mathobj.RS(Data[y][x].r.x,
                                               Data[y][x].r.y)[2])
    # convert to numpy array - just since matplotlib expects this format to plot
    RS_x_out = np.asarray(RS_x_data)
    RS_y_out = np.asarray(RS_y_data)
    RS_z_out = np.asarray(RS_z_data)
    # remove 0 elements from line plot (caused lines to jump to the origin when hit a singularity):
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            if RS_x_out[y][x] == 0 and RS_y_out[y][x] == 0 and RS_z_out[y][x] == 0 and x != 0 and y != 0:
                RS_x_out[y][x] = np.nan
                RS_y_out[y][x] = np.nan
                RS_z_out[y][x] = np.nan

    # now figure out how to export R's as json data
    with open('rx.json', 'w') as f:
        json.dump(R_x_data, f)
    with open('ry.json', 'w') as f:
        json.dump(R_y_data, f)
    with open('rz.json', 'w') as f:
        json.dump(R_z_data, f)
    with open('rx2.json', 'w') as f:
        json.dump(R_x_data2, f)
    with open('ry2.json', 'w') as f:
        json.dump(R_y_data2, f)
    with open('rz2.json', 'w') as f:
        json.dump(R_z_data2, f)
    with open('rsx.json', 'w') as f:
        json.dump(RS_x_data, f)
    with open('rsy.json', 'w') as f:
        json.dump(RS_y_data, f)
    with open('rsz.json', 'w') as f:
        json.dump(RS_z_data, f)
    # seems to work well
    return


def json_out_xy(mathobj, Data):
    # test this to ensure works, prints out x,y to be fed into RM, RS
    # now do not conver to numpy, json does not recognize...
    # first decompose to an array of r values as in plotting
    R_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(PHISIZE):
        for x in range(THETASIZE):
            # we then store the x, y, and z coordinates of the final surface at each point:
            R_x_data[y][x] = float(Data[y][x].romega.x)

            R_y_data[y][x] = float(Data[y][x].romega.y)

    # now extract 3d initial curves (vs above is for final curves):
    RS_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    RS_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]

    for y in range(PHISIZE):
        for x in range(THETASIZE):
            RS_x_data[y][x] = float(Data[y][x].r.x)
            RS_y_data[y][x] = float(Data[y][x].r.y)

    # now figure out how to export R's as json data
    with open('2drx.json', 'w') as f:
        json.dump(R_x_data, f)
    with open('2dry.json', 'w') as f:
        json.dump(R_y_data, f)
    with open('2drsx.json', 'w') as f:
        json.dump(RS_x_data, f)
    with open('2drsy.json', 'w') as f:
        json.dump(RS_y_data, f)
    # seems to work well
    return


def csv_xy(mathobj, Data):
    R_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(0, PHISIZE, int(PHISIZE/5)):
        for x in range(THETASIZE):
            R_x_data[y][x] = float(mathobj.R(Data[y][x].r[0].x, Data[y][x].r[0].y, 0)[
                0])
            R_y_data[y][x] = float(mathobj.R(Data[y][x].r[0].x,
                                             Data[y][x].r[0].y, 0)[1])
            file_name = "row"+str(y)+".csv"
            file = open(file_name, "a")
            file.write(str(R_x_data[y][x]) + "," +
                       str(R_y_data[y][x]) + ", 0\n")


def epsilon_csv_xy(Data):
    # prints for initial flat sheet with R being x,y,0... want to get epsilon lines too...
    R_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    R_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    # for epsilon thickness:
    eR_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    eR_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(0, PHISIZE, int(PHISIZE/5)):
        # clear out old data:
        file_name = "row"+str(y)+".csv"
        file = open(file_name, "w")
        file.write("")
        efile_name = "erow"+str(y)+".csv"
        efile = open(efile_name, "w")
        efile.write("")
        for x in range(THETASIZE):
            R_x_data[y][x] = float(Data[y][x].r[0].x)
            R_y_data[y][x] = float(Data[y][x].r[0].y)
            # want epsilon to be smaller than line spacing... how to do? maybe leave out lines until far enough away
            eR_x_data[y][x] = 0.5*STEP * \
                float(Data[y][x].n[0].perp().times(
                    Data[y][x].epsilon).x) + R_x_data[y][x]
            eR_y_data[y][x] = 0.5*STEP * \
                float(Data[y][x].n[0].perp().times(
                    Data[y][x].epsilon).y) + R_y_data[y][x]
            file_name = "row"+str(y)+".csv"
            file = open(file_name, "a")
            file.write(str(R_x_data[y][x]) + "," + str(R_y_data[y][x]) + "\n")
            efile_name = "erow"+str(y)+".csv"
            efile = open(efile_name, "a")
            efile.write(str(eR_x_data[y][x]) + "," +
                        str(eR_y_data[y][x]) + "\n")


def spaced_epsilon_csv_xy(Data, quadrant):
    # add variable line numbers later: jump_size = int(PHISIZE / lines)
    # now properly spaces cavity and bulk
    # prints for initial flat sheet with R being x,y,0... want to get epsilon lines too...
    #R_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    #R_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    # for epsilon thickness:
    #eR_x_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    #eR_y_data = [[0.0 for x in range(THETASIZE)] for y in range(PHISIZE)]
    for y in range(0, PHISIZE - int(PHISIZE/5), int(PHISIZE/5)):
        # ensure not at end ^
        # clear out old data:
        file_name = "q"+str(quadrant)+"row"+str(y)+".csv"
        file = open(file_name, "w")
        file.write("")
        efile_name = "q"+str(quadrant)+"erow"+str(y)+".csv"
        efile = open(efile_name, "w")
        efile.write("")
        for x in range(THETASIZE - 1):
            R_x = float(Data[y][x].r[0].x)
            R_y = float(Data[y][x].r[0].y)
            # difference between this and next row (do eg sqrt dx^2dy^2?):
            # d + dw
            # integrate beta dv, from vi to vf (next line)
            # trap. rule, for any grid spacing (though dv const-uniform, could use faster but ok)
            # get sum f(x) + f(f+1)/2 * dv
            # int_beta = 0
            # # might be bad since STEP maybe not true dv?
            # for step in range(y, y + int(PHISIZE/5) - 1, 1):
            #     int_beta += STEP*(Data[y + 1][x].beta +
            #                       Data[y][x].beta)/2  # dv = STEP?
            # d = Data[y][x].epsilon*int_beta

            R_xp = float(Data[y + int(PHISIZE/5)][x].r[0].x)
            R_yp = float(Data[y + int(PHISIZE/5)][x].r[0].y)
            dist = sqrt((R_yp - R_y)**2 + (R_xp - R_x)**2)
            d = Data[y][x].epsilon * dist
            # again take second line in nperp step?
            er_x = float(Data[y][x].n[0].perp().times(d).x) + R_x
            er_y = float(Data[y][x].n[0].perp().times(d).y) + R_y
            file_name = "q"+str(quadrant)+"row"+str(y)+".csv"
            file = open(file_name, "a")
            file.write(str(R_x) + "," + str(R_y) + "\n")
            efile_name = "q"+str(quadrant)+"erow"+str(y)+".csv"
            efile = open(efile_name, "a")
            efile.write(str(er_x) + "," +
                        str(er_y) + "\n")


def balanced_epsilon_csv_xy(Data, quadrant, num_lines):
    # phisize/8 i.e. ~8 lines
    # will do half an epsilon step in each direction, so cavities
    # are *centered* on lines, not starting on lines
    v_sgn = 1 if (quadrant == 1 or quadrant == 2) else -1
    for y in range(0, PHISIZE - int(PHISIZE/num_lines), int(PHISIZE/num_lines)):
        # ensure not at end ^
        # clear out old data:
        file_name = "q"+str(quadrant)+"row"+str(y)+".csv"
        file = open(file_name, "w")
        file.write("")
        efile_name = "q"+str(quadrant)+"erow"+str(y)+".csv"
        efile = open(efile_name, "w")
        efile.write("")
        e2file_name = "q"+str(quadrant)+"e2row"+str(y)+".csv"
        e2file = open(e2file_name, "w")
        e2file.write("")
        for x in range(THETASIZE - 1):
            R_x = float(Data[y][x].r[0].x)
            R_y = float(Data[y][x].r[0].y)

            R_xp = float(Data[y + int(PHISIZE/num_lines)][x].r[0].x)
            R_yp = float(Data[y + int(PHISIZE/num_lines)][x].r[0].y)
            dist = sqrt((R_yp - R_y)**2 + (R_xp - R_x)**2)
            d = Data[y][x].epsilon * dist * 0.5
            # again take second line in nperp step?
            er_x = float(Data[y][x].n[0].perp().times(v_sgn*d).x) + R_x
            er_y = float(Data[y][x].n[0].perp().times(v_sgn*d).y) + R_y
            if y >= int(PHISIZE/num_lines):
                # going opposite direction:
                R_xp2 = float(Data[y - int(PHISIZE/num_lines)][x].r[0].x)
                R_yp2 = float(Data[y - int(PHISIZE/num_lines)][x].r[0].y)
                dist2 = sqrt((R_yp2 - R_y)**2 + (R_xp2 - R_x)**2)
                d2 = Data[y][x].epsilon * dist2 * 0.5
                er_x2 = float(Data[y][x].n[0].perp().times(-v_sgn*d2).x) + R_x
                er_y2 = float(Data[y][x].n[0].perp().times(-v_sgn*d2).y) + R_y
                e2file_name = "q"+str(quadrant)+"e2row"+str(y)+".csv"
                e2file = open(e2file_name, "a")
                e2file.write(str(er_x2) + "," +
                             str(er_y2) + "\n")
            file_name = "q"+str(quadrant)+"row"+str(y)+".csv"
            file = open(file_name, "a")
            file.write(str(R_x) + "," + str(R_y) + "\n")
            efile_name = "q"+str(quadrant)+"erow"+str(y)+".csv"
            efile = open(efile_name, "a")
            efile.write(str(er_x) + "," +
                        str(er_y) + "\n")


def v_curve_balanced_epsilon_csv_xy(Data, quadrant, num_lines):
    # phisize/8 i.e. ~8 lines
    # will do half an epsilon step in each direction, so cavities
    # are *centered* on lines, not starting on lines
    v_sgn = 1 if (quadrant == 1 or quadrant == 2) else -1
    for x in range(0, THETASIZE - int(THETASIZE/num_lines), int(THETASIZE/num_lines)):
        # ensure not at end ^
        # clear out old data:
        file_name = "q"+str(quadrant)+"col"+str(x)+".csv"
        file = open(file_name, "w")
        file.write("")
        efile_name = "q"+str(quadrant)+"ecol"+str(x)+".csv"
        efile = open(efile_name, "w")
        efile.write("")
        e2file_name = "q"+str(quadrant)+"e2col"+str(x)+".csv"
        e2file = open(e2file_name, "w")
        e2file.write("")
        for y in range(PHISIZE - 1):
            R_x = float(Data[y][x].r[0].x)
            R_y = float(Data[y][x].r[0].y)

            R_xp = float(Data[y][x + int(THETASIZE/num_lines)].r[0].x)
            R_yp = float(Data[y][x + int(THETASIZE/num_lines)].r[0].y)
            dist = sqrt((R_yp - R_y)**2 + (R_xp - R_x)**2)
            d = Data[y][x].epsilon * dist * 0.5
            # again take second line in nperp step?
            er_x = float(Data[y][x].n[0].times(v_sgn*d).x) + R_x
            er_y = float(Data[y][x].n[0].times(v_sgn*d).y) + R_y
            if x >= int(THETASIZE/num_lines):
                # going opposite direction:
                R_xp2 = float(Data[y][x - int(THETASIZE/num_lines)].r[0].x)
                R_yp2 = float(Data[y][x - int(THETASIZE/num_lines)].r[0].y)
                dist2 = sqrt((R_yp2 - R_y)**2 + (R_xp2 - R_x)**2)
                d2 = Data[y][x].epsilon * dist2 * 0.5
                er_x2 = float(Data[y][x].n[0].times(-v_sgn*d2).x) + R_x
                er_y2 = float(Data[y][x].n[0].times(-v_sgn*d2).y) + R_y
                e2file_name = "q"+str(quadrant)+"e2col"+str(x)+".csv"
                e2file = open(e2file_name, "a")
                if ((x == 0 and y == 0) or not(er_x2 == 0 and er_y2 == 0)):
                    e2file.write(str(er_y2) + "," +
                                 str(er_x2) + "\n")
            file_name = "q"+str(quadrant)+"col"+str(x)+".csv"
            file = open(file_name, "a")
            if ((x == 0 and y == 0) or not (R_y == 0 and R_x == 0)):
                file.write(str(R_y) + "," + str(R_x) + "\n")
            efile_name = "q"+str(quadrant)+"ecol"+str(x)+".csv"
            efile = open(efile_name, "a")
            if ((x == 0 and y == 0) or not(er_y == 0 and er_x == 0)):
                efile.write(str(er_y) + "," +
                            str(er_x) + "\n")


def diagnostic_epsilon_csv_xy(Data, quadrant):
    # prints out all lines and just epsilon not offset
    file_name = "diagnostic_q"+str(quadrant)+".csv"
    file = open(file_name, "w")
    file.write("")  # clear old data
    file = open(file_name, "a")  # set to append for inserting data
    file.write(
        "u_x, u_y, r_x, r_y, n_x, n_y, alpha, beta, epsilon, b, s, p, q, d_u p\n")
    for y in range(0, PHISIZE):
        for x in range(THETASIZE):
            R_x = float(Data[y][x].r[0].x)
            R_y = float(Data[y][x].r[0].y)
            # now include epsilon, u, v too
            file.write(str(Data[y][x].uv.x) + ", " + str(Data[y][x].uv.y) + ", " + str(R_x) + "," + str(R_y) + ", " +
                       str(Data[y][x].n[0].x) + ", " + str(Data[y][x].n[0].y) + ", " +
                       str(Data[y][x].alpha) + "," + str(Data[y][x].beta) +
                       "," + str(Data[y][x].epsilon) + "," + str(Data[y][x].b) + "," +
                       str(Data[y][x].s) + "," + str(Data[y][x].p) + "," + str(Data[y][x].q) + "," + str(Data[y][x].dup) + "\n")
