#!/bin/bash

# Variables
LOCATION="CONUS"
BUCKET="noaa-mrms-pds"
REGION="us-east-1"
RADAR="RadarOnly_QPE_24H_00.00"
PASS1="MultiSensor_QPE_24H_Pass1_00.00"
PASS2="MultiSensor_QPE_24H_Pass2_00.00"
START_DATE=$(date -u -d "-4 years" +%Y-%m-%d)
END_DATE=$(date -u +%Y-%m-%d)

echo "Syncing MRMS QPE data for $LOCATION from $START_DATE to $END_DATE..."

# Loop from start date to end date
CURRENT_DATE="$START_DATE"
while [[ "$CURRENT_DATE" < "$END_DATE" || "$CURRENT_DATE" == "$END_DATE" ]]; do
  # Extract parts
  DDATE=$(date -u -d "$CURRENT_DATE" +%Y%m%d)
  MONTH=$(date -u -d "$CURRENT_DATE" +%m)
  DAY=$(date -u -d "$CURRENT_DATE" +%d)

  # S3 path and local destination
  SRC_PATH="s3://${BUCKET}/${LOCATION}/${RADAR}/${DDATE}/"
  DEST_DIR="/nfs/gce/projects/crocus/data/mrms/radar_qpe/${DDATE}/"

  PASS1_PATH="s3://${BUCKET}/${LOCATION}/${PASS1}/${DDATE}/"
  PASS1_DEST="/nfs/gce/projects/crocus/data/mrms/multisensor_qpe_pass1/${DDATE}/"

  PASS2_PATH="s3://${BUCKET}/${LOCATION}/${PASS2}/${DDATE}/"
  PASS2_DEST="/nfs/gce/projects/crocus/data/mrms/multisensor_qpe_pass2/${DDATE}/"

  # Create local directory if it doesn't exist
  mkdir -p "$DEST_DIR"

  echo "Syncing $CURRENT_DATE..."
  aws s3 sync "$SRC_PATH" "$DEST_DIR" \
    --region $REGION \
    --no-sign-request

  aws s3 sync "$PASS2_PATH" "$PASS2_DEST" \
    --region $REGION \
    --no-sign-request

  aws s3 sync "$PASS1_PATH" "$PASS1_DEST" \
    --region $REGION \
    --no-sign-request

  # Move to next day
  CURRENT_DATE=$(date -u -d "$CURRENT_DATE +1 day" +%Y-%m-%d)
done

echo "✅ All done syncing data for $LOCATION from the last 4 years."