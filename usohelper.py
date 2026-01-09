"""
Module and script with utilities for Thunderbrid South.



todo: observability chart for moon and target, with time in Bonn time

"""

from pathlib import Path
import argparse
import logging

from zoneinfo import ZoneInfo
from datetime import datetime, timedelta

import numpy as np


import astropy
import astropy.time
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, Angle
from pytz import timezone
from astroplan import Observer, FixedTarget


# Set up logging, adjust format
#logging_format = "%(levelname)s - %(message)s"
logging_format = "| %(message)s"
logging.basicConfig(format=logging_format, level=logging.INFO)
logger = logging.getLogger(__name__)


tzchile = ZoneInfo("America/Santiago")
tzvancouver = ZoneInfo("America/Vancouver")
tzbonn = ZoneInfo("Europe/Berlin")
#tzutc = ZoneInfo("UTC")


def format_dt(dt):
    """Format datetime object for display."""
    return dt.strftime("%a %Y-%m-%d %H:%M")


class AstroPlanWrapper():
    def __init__(self):

        # Define the observer location for Thunderbird South
        longitude = -70.8531 * u.deg
        latitude = -30.5261 * u.deg
        elevation = 1710 * u.m
        location = EarthLocation.from_geodetic(longitude, latitude, elevation)

        self.observer = Observer(name='Thunderbird South',
                                 location=location,
                                 pressure=1.0 * u.bar,
                                 relative_humidity=0.2,
                                 temperature=10 * u.deg_C,
                                 timezone=timezone('America/Santiago'),
                                 description="Thunderbird South")

        self.rng = np.random.default_rng()


    def night_overview(self, dt):
        """Compute and log an overview of sunset/rise and astro/nautical twilights for the night following the given datetime.
        
        dt is a datetime object.
        """

        def format_with_tz(time, tz):
            dt = self.observer.astropy_time_to_datetime(time)
            dt = dt.astimezone(tz)
            return dt.strftime("%a %Y-%m-%d %H:%M") + f" [{str(dt.tzinfo).split('/')[-1]}]"
        

        time = self.observer.datetime_to_astropy_time(dt)
        logger.info(f"=== Night of {format_with_tz(time, tzchile)} ===")
        
        sunset = self.observer.sun_set_time(time, which='next')
        nautstart = self.observer.twilight_evening_nautical(time, which='next')
        astrstart = self.observer.twilight_evening_astronomical(time, which='next')
        astrend = self.observer.twilight_morning_astronomical(time, which='next')
        nautend = self.observer.twilight_morning_nautical(time, which='next')
        sunrise = self.observer.sun_rise_time(time, which='next')

        logger.info(f"{'Sunset':<16}: {format_with_tz(sunset, tzchile)} = {format_with_tz(sunset, tzbonn)} = {format_with_tz(sunset, tzvancouver)}")
        logger.info(f"{'Nautical Start':<16}: {format_with_tz(nautstart, tzchile)} = {format_with_tz(nautstart, tzbonn)} = {format_with_tz(nautstart, tzvancouver)}")
        logger.info(f"{'Astron. Start':<16}: {format_with_tz(astrstart, tzchile)} = {format_with_tz(astrstart, tzbonn)} = {format_with_tz(astrstart, tzvancouver)}")
        logger.info(f"{'Astron. End':<16}: {format_with_tz(astrend, tzchile)} = {format_with_tz(astrend, tzbonn)} = {format_with_tz(astrend, tzvancouver)}")
        logger.info(f"{'Nautical End':<16}: {format_with_tz(nautend, tzchile)} = {format_with_tz(nautend, tzbonn)} = {format_with_tz(nautend, tzvancouver)}")
        logger.info(f"{'Sunrise':<16}: {format_with_tz(sunrise, tzchile)} = {format_with_tz(sunrise, tzbonn)} = {format_with_tz(sunrise, tzvancouver)}")


    def prepare_program(self, targetname, ra=None, dec=None, filter='R', exptime=60, nexp=5, stare=False, outputfile=None):
        """
        Prepares a "program" file with dithered exposures of the given target.
        
        :param targetname: name of the target
        :param ra: optional RA of target
        :param dec: optional Dec of target

        Reference output format:
        target: IC1613, ra: 01 04 48.00, dec: 02 07 04.0, filter: R, exposure time: 60, exposure type: light, rotator: absolute, image prefix: IC1613_R, focuser: 21200, nexposure: 1

        """
        dither = not stare
        dither_radius = 15 * u.arcsec

        if ra is not None and dec is not None:
            coordinates = SkyCoord(Angle(ra), Angle(dec), frame='icrs')
            target = FixedTarget(name=targetname, coord=coordinates)
        else: # We use Simbad to resolve the name
            target = FixedTarget.from_name(targetname)

        logger.info(f"Preparing program for target '{targetname}'")
        #logger.info(f"{str(target)}")

        targetname_safe = targetname.replace(" ", "_")

        output_lines = []
        output_lines.append(f"# {nexp} exposures of '{targetname}' in filter {filter} with {exptime} s each")
        for i in range(nexp):
            
            if dither:
                offset_along_ra = self.rng.uniform(low=-1.0, high=1.0) * dither_radius
                offset_along_dec = self.rng.uniform(low=-1.0, high=1.0) * dither_radius
                this_exp_coord = target.coord.spherical_offsets_by(offset_along_ra, offset_along_dec)
            else:
                this_exp_coord = target.coord

            line = f"target: {targetname_safe}, ra: {this_exp_coord.ra.to_string(unit=u.hour, sep=' ', pad=True, precision=2)}, dec: {this_exp_coord.dec.to_string(unit=u.degree, sep=' ', pad=True, precision=2, alwayssign=True)}, filter: {filter}, exposure time: {exptime}, exposure type: light, rotator: absolute, image prefix: {targetname_safe}_{filter}, focuser: 21200, nexposure: 1"

            output_lines.append(line)
        
        if outputfile is not None:
            with open(outputfile, 'w') as f:
                for line in output_lines:
                    f.write(line + "\n")
            logger.info(f"Wrote program to file: {outputfile}")
        else:
            print("\n".join(output_lines))



def current_time():
    """Return the current time in various time zones."""
    now_utc = datetime.now(tz=ZoneInfo("UTC"))
    output = {
        "UTC": now_utc,
        "Chile": now_utc.astimezone(tzchile),
        "Vancouver": now_utc.astimezone(tzvancouver),
        "Bonn": now_utc.astimezone(tzbonn),
    }
    
    logger.info("=== Current time ===")
    for zone, dt in output.items():
        logger.info(f"{zone:<12}: {format_dt(dt)} [{dt.tzinfo}]")

    return output



if __name__ == "__main__":
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Utilities for Thunderbrid South")
    parser.add_argument("-d", "--date", type=str, help="set a local date preceeding the night, in YYYY-MM-DD (default: today)", default=None)
    parser.add_argument("-v", "--nightoverview", action="store_true", help="print night overview")
    parser.add_argument("-t", "--targetname", type=str, help="name of target", default=None)
    parser.add_argument("--ra", type=str, help="RA of target (optional, in any format that can be parsed as Angle, e.g. '2h43m52s')", default=None)
    parser.add_argument("--dec", type=str, help="Dec of target (optional, in any format that can be parsed as Angle, e.g. '+25d30m')", default=None)
    parser.add_argument("-f", "--filter", type=str, help="filter (typically R, V, or B)", default='R')
    parser.add_argument("-e", "--exptime", type=int, help="exposure time in seconds (default: 60)", default=60)
    parser.add_argument("-n", "--nexp", type=int, help="number of exposures (default: 5)", default=5)
    parser.add_argument("-s", "--stare", action="store_true", help="'stare', i.e., do NOT dither the target")
    parser.add_argument("-p", "--prog", action="store_true", help="create program")
    parser.add_argument("-o", "--outputfile", type=str, help="output file for program", default=None)
    args = parser.parse_args()

    # Defining a local noon time at the observatory for the given date.
    # The "night" is the evening after that noon.
    if args.date is None:
        args.date = datetime.now(tz=tzchile).strftime("%Y-%m-%d")
    localnoon = datetime.strptime(args.date, "%Y-%m-%d")
    localnoon = datetime(localnoon.year, localnoon.month, localnoon.day, 12, 0, 0, tzinfo=tzchile)
    
    

    if args.nightoverview:
        ap = AstroPlanWrapper()
        ap.night_overview(localnoon)


    if args.prog:
        if args.targetname is None:
            logger.error("Please provide at least a target name.")
        else:
            ap = AstroPlanWrapper()
            ap.prepare_program(args.targetname, ra=args.ra, dec=args.dec, filter=args.filter, exptime=args.exptime, nexp=args.nexp, stare=args.stare, outputfile=args.outputfile)

    #current_time()
   

