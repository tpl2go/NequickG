""" Common sphereical trigonometry functions"""

import numpy as np

def greatcircle(lat1, lon1, lat2, lon2):
    """angle spanned by great arc connecting point 1 and 2
    symbol in reference document: delta"""
    DR = np.pi / 180.0

    # Spherical law of cosines
    cosine_delta = np.sin(lat1 * DR) * np.sin(lat2 * DR) + np.cos(lat1 * DR) * np.cos(lat2 * DR) * np.cos(
        (lon2 - lon1) * DR)
    sine_delta = np.sqrt(1 - cosine_delta ** 2)

    return np.arctan2(sine_delta, cosine_delta) * 180 / np.pi


def zenith(h1, lat1, lon1, h2, lat2, lon2):
        """
        zenith angle of point 2 from point 1
        :param lat1: [deg]
        :param lon1: [deg]
        :param lat2: [deg]
        :param lon2: [deg]
        :return: [deg]
        """
        delta = greatcircle(lat1, lon1, lat2, lon2)

        DR = np.pi / 180
        cosine = np.cos(delta * DR)
        sine = np.sin(delta * DR)

        r1 = 6371.2 + h1
        r2 = 6371.2 + h2

        zenith = np.arctan2(sine, cosine - r1 / r2) * 180.0 / np.pi

        return zenith

def azimuth(lat1, lon1, lat2, lon2):
        """ Imagine a bearing specified by a great arc from point 1 to point 2 on a
        unit sphere. Azimuth angle is the angle from the north bearing to current bearing"""

        delta = greatcircle(lat1, lon1, lat2, lon2)
        DR = np.pi / 180

        # spherical cosine rule on spherical triangle between pole and great circle
        cosine_sigma = (np.sin(lat2 * DR) - np.sin(lat1 * DR) * np.cos(delta * DR)) / (
            np.cos(lat1 * DR) * np.sin(delta * DR))
        # spherical sine rule on spherical triangle between pole and great circle
        sine_sigma = (np.cos(lat2 * DR) * np.sin((lon2 - lon1) * DR)) / np.sin(delta * DR)

        return np.arctan2(sine_sigma, cosine_sigma) * 180/np.pi
