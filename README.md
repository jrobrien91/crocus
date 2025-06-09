# crocus
Repo to store personal analysis for the Community Research on Climate and Urban Science ([CROCUS](https://crocus-urban.org/))

## To Do
| Topic | Description | Keywords | Which Paper/Poster? | Date Completed |
| :---: | :---------: | :------: | :-----------------: | :------------: |
| MRMS - Three Panal Plot | Specific Storm Analysis | MRMS; WXT | Precipitation Network Paper | [23 April 2025](https://github.com/jrobrien91/crocus/pull/20) |
| MRMS - 24 hr Precip Accum | Specific Storm Analysis | MRMS; WXT | Preciptiation Network Paper | [23 April 2025](https://github.com/jrobrien91/crocus/pull/20) |
| MRMS - 1 Hr Precip Accum | CROCUS Urban Flooding | MRMS; WXT | Preciptiation Network Paper | [8 May 2025](https://github.com/jrobrien91/crocus/pull/23) |
| MRMS - Direct Micronet Comparison | CROCUS Urban Flooding | MRMS; WXT | Preciptiation Network Paper | TBD |
| US Water API | USGS Rain Gauge for Precipitaiton Comparison | MRMS, WXT | Preciptiation Paper | TBD |
| CoCoRaHS API / Download | Need Gridded version for comparison | MRMS, WXT | Precipitaiton Paper | TBD |
| PBL Heights - Urban Canyons Radiosonde | Direct Radiosonde vs Lidar Comparison | Sonde, CL-61, DL | PBL Height Paper | TBD | 
| AQT 530 vs 560 Analysis | Performance of the AQT / Standars | AQT, QC | CROCUS - AEROMMA Paper | TBD |
| TEMPO N02 Analysis | Is the AQT capable of representing this? | AQT; QC | CROCUS - AEROMMA Paper | TBD |
| PANDORA Network API | Needed for Ozone over the city | AQ; AEROMMA | CROCUS - AEROMMA Paper | TBD |
| WXT Winter Analysis | Performance of the WXT in Snow/Mixed Precip | WXT; WinterWX | CROCUS Winter Precip | TBD |


## Necessary Workflow
### Switching between local and GCE
1. Track a Remote Github Branch
```bash
git branch --track branch-name origin/branch-name
```
### Updating cronjobs on GCE
1. To update the crontab
```bash
crontab -e
```

2. How to backup crontab
```bash
crontab -l > gce_crontab.bak
```
### Transfering data between locations using tmux and rsync
1. Start a new tmux session
   ```bash
   tmux new -s {session_name}
   ```
2. Within the tmux session, copy data
   ```bash
   rsync -avz user@computer:path_to_data path_to_save_location
   ```
3. Detach from the tmux session
   ```bash
   CTRL+B, D
   ```
4. Reattach to the tmux session when needed
   ```bash
   tmux attach -t {session_name}
   ```
5. Kill a tmux session
   ```bash
   tmux kill-session -t {session_name}
   ```
