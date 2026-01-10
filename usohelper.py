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

import matplotlib.pyplot as plt

import astropy
import astropy.time
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, Angle
from pytz import timezone
from astroplan import Observer, FixedTarget
from astroplan.plots import plot_airmass


# Set up logging, adjust format
#logging_format = "%(levelname)s - %(message)s"
logging_format = "| %(message)s"
logging.basicConfig(format=logging_format, level=logging.INFO)
logger = logging.getLogger(__name__)


tzchile = ZoneInfo("America/Santiago")
tzvancouver = ZoneInfo("America/Vancouver")
tzbonn = ZoneInfo("Europe/Berlin")
tzutc = ZoneInfo("UTC")

custom_tz_names = {
        'America/Santiago': 'Chile',
        'America/Vancouver': 'Vancouver',
        'Europe/Berlin': 'Bonn',
        'UTC': 'UTC',
    }


def format_dt(dt, long=True):
    """Format datetime object for display."""
    if long:
        return dt.strftime("%a %Y-%m-%d %H:%M")
    else:
        return dt.strftime("%a %d %H:%M")


def format_dt_with_tz(dt, tz, **kwargs):
    dt_with_tz = dt.astimezone(tz)
    tz_name = str(dt_with_tz.tzinfo)
    return format_dt(dt_with_tz, **kwargs) + f" [{custom_tz_names.get(tz_name, tz_name)}]"

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

        self.rng = np.random.default_rng(seed=42)


    def night_overview(self, dt):
        """Compute and log an overview of sunset/rise and astro/nautical twilights for the night following the given datetime.
        
        dt is a datetime object.
        """

        def format(time, tz, **kwargs):
            dt = self.observer.astropy_time_to_datetime(time)
            return format_dt_with_tz(dt, tz, **kwargs)
            
        
        time = self.observer.datetime_to_astropy_time(dt)
        logger.info(f"=== Night of {format(time, tzchile)} ===")
        
        sunset = self.observer.sun_set_time(time, which='next')
        nautstart = self.observer.twilight_evening_nautical(time, which='next')
        astrstart = self.observer.twilight_evening_astronomical(time, which='next')
        astrend = self.observer.twilight_morning_astronomical(time, which='next')
        nautend = self.observer.twilight_morning_nautical(time, which='next')
        sunrise = self.observer.sun_rise_time(time, which='next')

        logger.info(f"{'Sunset':<16}: {format(sunset, tzchile, long=False)} = {format(sunset, tzbonn, long=False)} = {format(sunset, tzvancouver, long=False)} = {format(sunset, tzutc, long=False)}")
        logger.info(f"{'Nautical Start':<16}: {format(nautstart, tzchile, long=False)} = {format(nautstart, tzbonn, long=False)} = {format(nautstart, tzvancouver, long=False)} = {format(nautstart, tzutc, long=False)}")
        logger.info(f"{'Astron. Start':<16}: {format(astrstart, tzchile, long=False)} = {format(astrstart, tzbonn, long=False)} = {format(astrstart, tzvancouver, long=False)} = {format(astrstart, tzutc, long=False)}")
        logger.info(f"{'Astron. End':<16}: {format(astrend, tzchile, long=False)} = {format(astrend, tzbonn, long=False)} = {format(astrend, tzvancouver, long=False)} = {format(astrend, tzutc, long=False)}")
        logger.info(f"{'Nautical End':<16}: {format(nautend, tzchile, long=False)} = {format(nautend, tzbonn, long=False)} = {format(nautend, tzvancouver, long=False)} = {format(nautend, tzutc, long=False)}")
        logger.info(f"{'Sunrise':<16}: {format(sunrise, tzchile, long=False)} = {format(sunrise, tzbonn, long=False)} = {format(sunrise, tzvancouver, long=False)} = {format(sunrise, tzutc, long=False)}")


    def set_target(self, targetname, ra=None, dec=None, coords=None):

        if targetname is None:
            logger.error("Always provide at least a target name!")

        if ra is None and dec is None and coords is None:
            # Then we resolve the name via Simbad
            self.target = FixedTarget.from_name(targetname)

        elif coords is not None:
            # We parse the given coord string
            #ra_str, dec_str = coord.split()
            #coordinates = SkyCoord(Angle(ra_str), Angle(dec_str), frame='icrs')
            coordinates = SkyCoord(coords, frame='icrs')
            self.target = FixedTarget(name=targetname, coord=coordinates)

        elif ra is not None and dec is not None:
            coordinates = SkyCoord(Angle(ra), Angle(dec), frame='icrs')
            self.target = FixedTarget(name=targetname, coord=coordinates)
        
        else:
            logger.error("Either provide coords, or both ra and dec, or provide neither to resolve the name via Simbad!")

        self.targetname_safe = targetname.replace(" ", "_")
        
        

    def prepare_program(self, filter='R', exptime=60, nexp=5, stare=False, outputfile=None):
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

        logger.info(f"Preparing program for target '{self.target.name}'")
        
        output_lines = []
        output_lines.append(f"# {nexp} exposures of '{self.target.name}' in filter {filter} with {exptime} s each")
        for i in range(nexp):
            
            if dither:
                offset_along_ra = self.rng.uniform(low=-1.0, high=1.0) * dither_radius
                offset_along_dec = self.rng.uniform(low=-1.0, high=1.0) * dither_radius
                this_exp_coord = self.target.coord.spherical_offsets_by(offset_along_ra, offset_along_dec)
            else:
                this_exp_coord = self.target.coord

            line = f"target: {self.targetname_safe}, ra: {this_exp_coord.ra.to_string(unit=u.hour, sep=' ', pad=True, precision=2)}, dec: {this_exp_coord.dec.to_string(unit=u.degree, sep=' ', pad=True, precision=2, alwayssign=True)}, filter: {filter}, exposure time: {exptime}, exposure type: light, rotator: absolute, image prefix: {self.targetname_safe}_{filter}, focuser: 21200, nexposure: 1"

            output_lines.append(line)
        
        if outputfile is not None:
            with open(outputfile, 'w') as f:
                for line in output_lines:
                    f.write(line + "\n")
            logger.info(f"Wrote program to file: {outputfile}")
        else:
            print("\n".join(output_lines))


    def plot_observability_chart(self, dt):
        """
        """
        time = self.observer.datetime_to_astropy_time(dt)

        plot_airmass(self.target, self.observer, time, brightness_shading=True)
        plt.title(f"Observability of {self.target.name}")
        plt.tight_layout
        plt.show()


def current_time():
    """Return the current time in various time zones."""
    now = datetime.now(tz=tzutc)
    logger.info("=== Current time ===")
    logger.info(f"{format_dt_with_tz(now, tzchile, long=False)} = {format_dt_with_tz(now, tzbonn, long=False)} = {format_dt_with_tz(now, tzvancouver, long=False)} = {format_dt_with_tz(now, tzutc, long=False)}\n")
    



if __name__ == "__main__":
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Utilities for Thunderbrid South")
    parser.add_argument("-d", "--date", type=str, help="set a local date preceeding the night, in YYYY-MM-DD (default: today)", default=None)
    parser.add_argument("-v", "--nightoverview", action="store_true", help="print night overview")
    parser.add_argument("-t", "--targetname", type=str, help="name of target", default=None)
    parser.add_argument("--ra", type=str, help="RA of target (optional, in any format that can be parsed as Angle, e.g. '2h43m52s')", default=None)
    parser.add_argument("--dec", type=str, help="Dec of target (optional, in any format that can be parsed as Angle, e.g. '+25d30m')", default=None)
    parser.add_argument("--coords", type=str, help="Coordinates of target, in format 'RA DEC' (e.g. '2h43m52s +25d30m')", default=None)
    parser.add_argument("-f", "--filter", type=str, help="filter (typically R, V, or B)", default='R')
    parser.add_argument("-e", "--exptime", type=int, help="exposure time in seconds (default: 60)", default=60)
    parser.add_argument("-n", "--nexp", type=int, help="number of exposures (default: 5)", default=5)
    parser.add_argument("-s", "--stare", action="store_true", help="'stare', i.e., do NOT dither the target")
    parser.add_argument("-p", "--prog", action="store_true", help="create program")
    parser.add_argument("-o", "--outputfile", type=str, help="output file for program", default=None)
    parser.add_argument("-a", "--plotairmass", action="store_true", help="show observability chart for target")
    args = parser.parse_args()

    # Defining a local noon time at the observatory for the given date.
    # The "night" is the evening after that noon.
    if args.date is None:
        args.date = datetime.now(tz=tzchile).strftime("%Y-%m-%d")
    localnoon = datetime.strptime(args.date, "%Y-%m-%d")
    localnoon = datetime(localnoon.year, localnoon.month, localnoon.day, 12, 0, 0, tzinfo=tzchile)
    localmidnight = datetime(localnoon.year, localnoon.month, localnoon.day, 23, 59, 59, tzinfo=tzchile)
    
    # We allways print the current time
    current_time()

    if args.nightoverview:
        ap = AstroPlanWrapper()
        ap.night_overview(localnoon)

    if args.prog:
        ap = AstroPlanWrapper()
        ap.set_target(args.targetname, ra=args.ra, dec=args.dec, coords=args.coords)
        ap.prepare_program(filter=args.filter, exptime=args.exptime, nexp=args.nexp, stare=args.stare, outputfile=args.outputfile)

    if args.plotairmass:
        ap = AstroPlanWrapper()
        ap.set_target(args.targetname, ra=args.ra, dec=args.dec, coords=args.coords)
        ap.plot_observability_chart(localmidnight)
    
   

