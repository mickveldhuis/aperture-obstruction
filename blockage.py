import argparse

from src.aperture import TelescopeAperture, GuiderAperture, FinderAperture

parser = argparse.ArgumentParser(
            allow_abbrev=True, 
            description='Compute the % obstruction of the main aperture/finder/guider by the dome'
        )

# Add arguments
parser.add_argument('--az', action='store', type=float, default=0.0, help='dome azimuth: 0 to 360 deg | default: 0 deg')
parser.add_argument('--ha', action='store', type=float, default=0.0, help='telescope hour angle: -180 to 180 deg | default: 0 deg')
parser.add_argument('--dec', action='store', type=float, default=0.0, help='telescope declination -90 to 90 deg | default: 0 deg')
parser.add_argument('-a', '--aperture', action='store', type=str, default='telescope', help='select aperture: telescope, finder, guider | default: telescope')
parser.add_argument('-r', '--rate', action='store', type=int, default=100, help='no. rays (for decent results >100; preferably 1000+) | default: 100')

args = parser.parse_args()

if __name__ == '__main__':    
    blockage = None

    if args.aperture == 'telescope':
        telescope = TelescopeAperture(rate=args.rate)
        blockage = telescope.obstruction(args.ha, args.dec, args.az)

    elif args.aperture == 'guider':
        guider = GuiderAperture(rate=args.rate)
        blockage = guider.obstruction(args.ha, args.dec, args.az)

    elif args.aperture == 'finder':
        finder = FinderAperture(rate=args.rate)
        blockage = finder.obstruction(args.ha, args.dec, args.az)

    if blockage is not None:
        print('obstruction = {:.2%}'.format(blockage))
    else:
        print('the % obstruction could not be computed..!')