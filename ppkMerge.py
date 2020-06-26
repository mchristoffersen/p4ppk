import pandas as pd
import numpy as np
import sys, argparse
import matplotlib.pyplot as plt
from pathlib import Path


def vmag(v0):
    return np.sqrt(v0[:, 0] ** 2 + v0[:, 1] ** 2 + v0[:, 2] ** 2)


def vunit(v0):
    mag = vmag(v0)
    v1 = v0[:]

    v1[:, 0] /= mag
    v1[:, 1] /= mag
    v1[:, 2] /= mag

    return v1


def nevVec(lat, lon):
    # Constructs north, east, and up vectors at
    # given locations
    # Lat and lon in degrees
    lat = np.radians(lat)
    lon = np.radians(lon)

    un = np.column_stack(
        (-np.sin(lat) * np.cos(lon), -np.sin(lat) * np.sin(lon), np.cos(lat))
    )
    ue = np.column_stack((-np.sin(lon), np.cos(lon), np.zeros_like(lon)))
    uv = np.cross(ue, un)

    return [un, ue, uv]


def lle2ecef(lat, lon, elev):
    # For WGS84 ellipsoid
    # Elev is elev above ellipsoid
    # Lat and lon in degrees
    lat = np.radians(lat)
    lon = np.radians(lon)

    a = 6378137  # WGS84 major axis
    f = 1 / 298.257223563  # WGS84 flattening
    b = a - f * a  # WGS84 minor axis
    e = np.sqrt(1 - np.power(b / a, 2))

    N = a / np.sqrt(1 - (np.power(e, 2) * np.power(np.sin(lat), 2)))

    x = (N + elev) * np.cos(lat) * np.cos(lon)
    y = (N + elev) * np.cos(lat) * np.sin(lon)
    z = ((np.power(b / a, 2) * N) + elev) * np.sin(lat)

    return [x, y, z]


def ecef2lle(x, y, z):
    # From Olsen (1996)
    # For WGS84 ellipsoid
    # X, Y, Z in meters
    # Lat and lon returned in degrees
    # Elev is height above wgs84 ellipsoid

    a = 6378137  # WGS84 major axis
    f = 1 / 298.257223563  # WGS84 flattening
    e2 = 2 * f - np.power(f, 2)  # eccentricity^2

    a1 = a * e2
    a2 = np.power(a, 2) * np.power(e2, 2)
    a3 = a * np.power(e2, 2) / 2
    a4 = (5.0 / 2) * np.power(a, 2) * np.power(e2, 2)
    a5 = a1 + a3
    a6 = 1 - e2

    zp = np.abs(z)
    w2 = np.power(x, 2) + np.power(y, 2)
    w = np.sqrt(w2)
    z2 = np.power(z, 2)
    r2 = w2 + z2
    r = np.sqrt(r2)

    lon = np.arctan2(y, x)
    s2 = z2 / r2
    c2 = w2 / r2
    u = a2 / r
    v = a3 - a4 / r

    s = np.zeros(len(x))
    lat = np.zeros(len(x))
    ss = np.zeros(len(x))
    c = np.zeros(len(s))

    lp3 = c2 > 0.3
    s[lp3] = (zp[lp3] / r[lp3]) * (
        1 + c2[lp3] * (a1 + u[lp3] + s2[lp3] * v[lp3]) / r[lp3]
    )
    lat[lp3] = np.arcsin(s[lp3])
    ss[lp3] = np.power(s[lp3], 2)
    c[lp3] = np.sqrt(1 - ss[lp3])

    gep3 = np.logical_not(lp3)
    c[gep3] = (w[gep3] / r[gep3]) * (
        1 - s2[gep3] * (a5 - u[gep3] - c2[gep3] * v[gep3]) / r[gep3]
    )
    lat[gep3] = np.arccos(c[gep3])
    ss[gep3] = 1 - np.power(c[gep3], 2)
    s[gep3] = np.sqrt(ss[gep3])

    g = 1 - e2 * ss
    rg = a / np.sqrt(g)
    rf = a6 * rg
    u = w - rg * c
    v = zp - rf * s
    f = c * u + s * v
    m = c * v - s * u
    p = m / ((rf / g) + f)
    lat = lat + p
    elev = f + m * p / 2

    lat[z < 0] = -1 * lat[z < 0]

    lat = np.degrees(lat)
    lon = np.degrees(lon)

    lat[r < 100000] = np.nan
    lon[r < 100000] = np.nan
    elev[r < 100000] = np.nan

    return [lat, lon, elev]


def cli():
    parser = argparse.ArgumentParser(
        description="Apply RTKLIB PPK solution to DJI Phantom 4 RTK Surveys"
    )
    parser.add_argument("ProjDir", help="Survey Directory")
    args = parser.parse_args()
    return Path(args.ProjDir)


def main():
    path = cli()
    proj = path.parts[-1]
    tstamp = path.joinpath(proj + "_Timestamp.MRK")
    fix = path.joinpath(proj + "_PPKRAW.pos")
    stamp = path.joinpath(proj + "_PPKStamp.csv")

    ## Read in the photo timestamps
    TScols = [
        "seq",
        "sow",
        "week",
        "ofstN",
        "ofstE",
        "ofstV",
        "lat",
        "lon",
        "elev",
        "stdN",
        "stdE",
        "stdV",
        "status",
    ]
    ts = pd.read_csv(tstamp, delim_whitespace=True, names=TScols)

    # Clean up some of the fields
    ts["lat"] = ts["lat"].str.rstrip(",Lat").astype(np.float64)
    ts["lon"] = ts["lon"].str.rstrip(",Lon").astype(np.float64)
    ts["elev"] = ts["elev"].str.rstrip(",Ellh").astype(np.float64)

    # Clean up and convert to meters
    ts["ofstN"] = ts["ofstN"].str.rstrip(",N").astype(np.int32) / 1000.0
    ts["ofstE"] = ts["ofstE"].str.rstrip(",E").astype(np.int32) / 1000.0
    ts["ofstV"] = ts["ofstV"].str.rstrip(",V").astype(np.int32) / 1000.0

    ## Read in the PPK fix
    FXcols = [
        "date",
        "time",
        "lat",
        "lon",
        "elev",
        "Q",
        "ns",
        "sdn",
        "sde",
        "sdu",
        "sdne",
        "sdeu",
        "sdun",
        "age",
        "ratio",
    ]
    fx = pd.read_csv(fix, delim_whitespace=True, names=FXcols, comment='%')

    fx["gpst"] = pd.to_datetime(fx["date"] + " " + fx["time"])

    #plt.scatter(ts["lat"], ts["lon"])

    # Convert the PPK solution times to seconds of week
    fx["sow"] = np.zeros(len(fx["gpst"])).astype(np.float64)
    for i in range(len(fx["gpst"])):
        wd = fx["gpst"][i].weekday()
        hr = fx["gpst"][i].hour
        mn = fx["gpst"][i].minute
        sc = fx["gpst"][i].second
        ms = fx["gpst"][i].microsecond
        fx.loc[i, "sow"] = ((wd + 1) % 7) * 86400 + hr * 3600 + mn * 60 + sc + ms * 1e-6

    # Interpolate new drone GPS antenna locations
    lat_ppk = np.interp(ts["sow"], fx["sow"], fx["lat"])
    lon_ppk = np.interp(ts["sow"], fx["sow"], fx["lon"])
    elev_ppk = np.interp(ts["sow"], fx["sow"], fx["elev"])

    # Calculate unit north, east, and up vectors for each position
    un, ue, uv = nevVec(lat_ppk, lon_ppk)

    # Apply camera offsets to PPK solutions
    x_ppk, y_ppk, z_ppk = lle2ecef(lat_ppk, lon_ppk, elev_ppk)
    cor = (
        un * ts["ofstN"][:, np.newaxis]
        + ue * ts["ofstE"][:, np.newaxis]
        + uv * ts["ofstV"][:, np.newaxis]
    )
    x_ppk += cor[:, 0]
    y_ppk += cor[:, 1]
    z_ppk += cor[:, 2]

    lat_ppk, lon_ppk, elev_ppk = ecef2lle(x_ppk, y_ppk, z_ppk)

    #plt.scatter(lat_ppk, lon_ppk)
    #plt.show()

    out = np.column_stack((ts["seq"], lat_ppk, lon_ppk, elev_ppk))
    hdr = "photo,latitude,longitude,elevation"
    fstr = proj + "_%04d,%.10f,%.10f,%.4f"
    np.savetxt(stamp, out, header=hdr, delimiter=",", fmt=fstr)

main()
