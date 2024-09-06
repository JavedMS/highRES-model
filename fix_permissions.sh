#!/bin/bash
#cd /cluster/work/projects/ec85/joint-wind/model-aggregated
#cd shared_input
echo "Script started"
chgrp -R ec85-member-group .
find . -type f -exec chmod 770 {} +
find . -type d -exec chmod 2770 {} +
