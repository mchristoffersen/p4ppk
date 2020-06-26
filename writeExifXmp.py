import pandas as pd
import piexif, argparse
from pathlib import Path
from libxmp import XMPFiles, consts
import sys, pprint

def cli():
    parser = argparse.ArgumentParser(
        description="Write PPK corrected coordinates to JPG images"
    )
    parser.add_argument("ProjDir", help="Survey Directory")
    args = parser.parse_args()
    return Path(args.ProjDir)


def main():
    path = cli()
    proj = path.parts[-1]
    stamp = path.joinpath(proj + "_PPKStamp.csv")

    cols = ["name","lat","lon","elev"]
    st = pd.read_csv(stamp, names=cols, comment='#')
    
    for i in range(len(st)):
        jpg = path.joinpath(st["name"][i] + ".JPG")
        exif = piexif.load(str(jpg))

        xmpfile = XMPFiles(file_path=str(jpg), open_forupdate=True)
        xmp = xmpfile.get_xmp()
        DJIns = xmp.get_namespace_for_prefix("drone-dji")

        xmp.set_property(DJIns, 'AbsoluteAltitude', "+%.4f" % st["elev"][i])
        xmp.set_property(DJIns, 'GpsLatitude', "%.10f" % st["lat"][i])
        xmp.set_property(DJIns, 'GpsLongitude', "%.10f" % st["lon"][i])
        #xmp.set_property(DJIns, 'GpsLongtitude', "+%.10f" % st["lon"][i])

        xmpfile.put_xmp(xmp)
        xmpfile.close_file()

        lat = abs(st["lat"][i])
        latd = int(lat)
        latm = int((lat-latd)*60)
        lats = (lat-latd-latm/60)*3600

        lon = abs(st["lon"][i])
        lond = int(lon)
        lonm = int((lon-lond)*60)
        lons = (lon-lond-lonm/60)*3600

        exif["GPS"][2] = ((latd, 1), (latm, 1), (int(lats*1000000), 1000000))
        exif["GPS"][4] = ((lond, 1), (lonm, 1), (int(lons*1000000), 1000000))
        exif["GPS"][6] = (int(st["elev"][i]*1000), 1000)

        exifBytes = piexif.dump(exif)
        piexif.insert(exifBytes, str(jpg))

main()

