import os
import re
import fnmatch
import glob
from trackhub import Track, default_hub
from trackhub.upload import upload_hub, upload_track
from logging import getLogger

logger = getLogger()


hub, genomes_file, genome, trackdb = default_hub(
    hub_name="hisat2",
    genome="dm6",
    short_label="Flybase example Hisat2",
    long_label="Flybase example Hisat2",
    email="justin.fear@nih.gov")

# publicly accessible hub URL
hub.url = "http://helix.nih.gov/~fearjm/ncbi_remap/hubs/flybase_example_hub.txt"

# hub's location on remote host, for use with rsync
hub.remote_fn = "/data/fearjm/datashare/ncbi_remap/hubs/flybase_example_hub.txt"

# Make tracks for all bigWigs in current dir
for root, dirnames, filenames in os.walk('../output/alignment/raw'):
    for fn in fnmatch.filter(filenames, '*.bw'):
        regex = re.compile('../output/alignment/raw/(?P<experiment>.*?)/(?P<sample>.*?)/(?P=sample).fq.bam.(?P<strand>first|second).bw')
        m = re.match(regex, os.path.join(root, fn))
        label = '{experiment}_{sample}_{strand}'.format(**m.groupdict())

        # Parameters are checked for valid values, see
        track = Track(
            name=label,
            short_label=label,
            long_label=label,
            autoScale='off',
            local_fn=os.path.join(root, fn),
            tracktype='bigWig',
            )
        trackdb.add_tracks(track)

# Demonstrate some post-creation adjustments...here, just make control
# samples gray
for track in trackdb.tracks:
    if 'first' in track.name:
        track.add_params(color="255,0,0")
    elif 'second' in track.name:
        track.add_params(color="0,0,255")

# Render the hub to text files
hub.render()

# Upload the hub files and all the bigwig files using rsync.
kwargs = dict(host='biowulf.nih.gov', user='fearjm')
upload_hub(hub=hub, run_local=True, **kwargs)
for track, level in hub.leaves(Track):
    upload_track(track=track, run_local=True, **kwargs)
