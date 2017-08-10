import csv
import CCIR_MoDIP.modip
import numpy as np
from aux import *
import matplotlib.pyplot as plt


def solar_declination():
    """

    Compute sin(delta_Sun ), cos(delta_Sun ), the sine and cosine of the solar declination.
    :param month: [mth]
    :param universal_time: [hours and decimals]
    :return:(Cosine, Sine)
    """

    month = np.arange(1,13)
    universal_time = 12

    # Compute day of year at the middle of the month
    dy = 30.5 * month - 15
    # Compute time [days]:
    t = dy + (18 - universal_time) / 24
    # Compute the argument
    a_m = (0.9856 * t - 3.289) * (np.pi / 180)  # radians
    a_l = a_m + (1.916 * np.sin(a_m) + 0.020 * np.cos(2 * a_m) + 282.634) * (np.pi / 180)  # radians

    # Compute sine and cosine of solar declination
    Sine = 0.39782 * np.sin(a_l)
    Cosine = np.sqrt(1 - Sine ** 2)

    return (Cosine, Sine)

(Cosine, Sine) = solar_declination()
plt.plot(Cosine)
plt.plot(Sine)
plt.show()

#
# def modip(latitude, longitude):
#
#     # stModip = self.__read_stMoDIP__()
#     stModip = CCIR_MoDIP.modip.stModip
#
#     lngp = 36
#     dlatp = 5
#     dlngp = 10
#
#     if (latitude > 90) or (latitude < -90):
#         mu = 90
#         return mu
#
#     lng1 = (longitude + 180.0) / dlngp
#     sj = int(lng1) - 2
#     dj= lng1 - int(lng1)
#
#     if (sj < 0):
#         sj = sj + lngp
#     if (sj > (lngp - 3)):
#         sj = sj - lngp
#
#     lat1 = (latitude + 90.0) / dlatp + 1
#     si = int(lat1 - 1e-6) - 2
#     di = lat1 - si - 2
#
#
#     z = []
#     for k in range(1,5):
#         z_row = []
#         for j in range(1,5):
#             z_jk = stModip[si + j][sj+k + 1]
#             z_row.append(z_jk)
#         z_row.append(di)
#         z.append(interpolate(*z_row))
#     z.append(dj)
#     return interpolate(*z)
#
#
# #
# # A = np.zeros([180,360])
# #
# # for i in range(360):
# #     for j in range(180):
# #         A[j,i] = modip(j-90,i-180)
# #
# #
# # B = A[0::5,0::10]
# # print B
# # print
# # print np.array(CCIR_MoDIP.modip.stModip_nowrap)
# # print
# # print  np.array(CCIR_MoDIP.modip.stModip_nowrap)[:-1,:-1] - B
# # plt.imshow(A)
# # plt.figure()
# # plt.imshow(CCIR_MoDIP.modip.stModip_nowrap)
# # plt.show()
# #
#
# def modip2(latitude, longitude):
#
#     lon_index = int((longitude + 180)/10.0)
#
#     lat_index = int((latitude + 90)/5.0)
#
#     lon_excess = (longitude + 180)/10.0 - lon_index
#     lat_excess = (latitude + 90)/5.0 - lat_index
#
#     # print lat_excess, lon_excess
#
#     stModipnw = CCIR_MoDIP.modip.stModip_nowrap
#
#     # interpolate vertically first then horizontally. it shouldnt matter right?
#
#     square = np.empty([4,4])
#     for i in range(4):
#         for j in range(4):
#             square[j,i] = stModipnw[lat_index -1 + j][lon_index -1 + i]
#
#     row = []
#     for i in range(4):
#         row.append(interpolate(square[i,0],square[i,1],square[i,2],square[i,3],lon_excess))
#         # row.append(interpolate(square[0,i],square[1,i],square[2,i],square[3,i],lat_excess))
#
#     row.append(lat_excess)
#     # row.append(lon_excess)
#     return interpolate(*row)
#
#

#
# # print Square(1,1)
# # print modip(1,1)
# # print modip2(1,1)
# for i in range(6):
#     for j in range(11):
#         print modip(i,j)
#         print modip2(i,j)
#         print
