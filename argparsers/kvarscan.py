import ruffus as R


def parse_args():
    parser = R.cmdline.get_argparse(
        description="konnector",
        usage='require python-2.7.x',
        version='0.1')

    parser.add_argument('--bams', nargs='+', required=True)
    options = parser.parse_args()
    return options

OPTIONS=parse_args()
