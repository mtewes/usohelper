"""
Module and script with utilities for Thunderbrid South.



todo: observability chart for moon and target, with time in Bonn time

"""

from pathlib import Path
import argparse
import logging

from zoneinfo import ZoneInfo
from datetime import datetime, timedelta


import astropy
import astropy.time
import astropy.units as u
from astropy.coordinates import EarthLocation
from pytz import timezone
from astroplan import Observer

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


def compute_night(date):
    """Compute nautical and astro night times for a given date."""
    dt = datetime(date.year, date.month, date.day, 12, 0, 0, tzinfo=tzchile)

    # Nautical twilight
    nautical_start = dt - timedelta(hours=6, minutes=30)
    nautical_end = dt + timedelta(hours=6, minutes=30)

    # Astronomical twilight
    astro_start = dt - timedelta(hours=7, minutes=30)
    astro_end = dt + timedelta(hours=7, minutes=30)

    return {
        "nautical": (nautical_start, nautical_end),
        "astronomical": (astro_start, astro_end),
    }

if __name__ == "__main__":
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Utilities for Thunderbrid South")
    parser.add_argument("-d", "--date", type=str, help="set a local date preceeding the night, in YYYY-MM-DD (default: today)", default=None)
    parser.add_argument("-n", "--n", action="store_true", help="print night overview")
    args = parser.parse_args()

    # Defining a local noon time at the observatory for the given date.
    # The "night" is the evening after that noon.
    if args.date is None:
        args.date = datetime.now(tz=tzchile).strftime("%Y-%m-%d")
    localnoon = datetime.strptime(args.date, "%Y-%m-%d")
    localnoon = datetime(localnoon.year, localnoon.month, localnoon.day, 12, 0, 0, tzinfo=tzchile)
    
    

    if args.n:
        ap = AstroPlanWrapper()
        ap.night_overview(localnoon)


    #current_time()
   

