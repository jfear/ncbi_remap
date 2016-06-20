#!/usr/bin/env python
""" Takes all BigWig files and make a trackhub. """
from trackhub import Track, default_hub
from trackhub.upload import upload_hub, upload_track


if __name__ == '__main__':
    hub, genomes_file, genome, trackdb = default_hub(
            hub_name='hisat2_318_golden_samples',
            genome='dm6',
            short_label='hisat2 318 golden samples',
            long_label='hisat2 318 golden samples',
            email='justin.fear@nih.gov')

    # publicly accessible hub URL
    hub.url = 'http://helix.nih.gov/~fearjm/ncbi_remap/hub.txt'

    # hub's location on remote host
    hub.remote_fn = '/data/fearjm/datashare/ncbi_remap/hub.txt'

    # Make tracks for all bigWigs
    import glob, os
    for fn in glob.glob('*.bw'):
        label = fn.replace('.bw', '')

        track = Track(
                name=label,
                short_label=label
                long_label=label
                autoScale='off',
                local_fn=fn,
                tracktype='bigWig',
                )
        trackdb.add_tracks(track)

    # Demonstrate some post-creation adjustment
    for track in trackdb.tracks:
        if 'control' in track.name:
            track.add_params(color="100,100,100")

    hub.render()

    # Upload the hub files and all the bigwig files using rsync.
    kwargs = dict(host='XXX', user='fearjm')
    upload_hub(hub=hub, **kwargs)
    for track, level in hub.leaves(Track):
        upload_track(track=track, **kwargs)
