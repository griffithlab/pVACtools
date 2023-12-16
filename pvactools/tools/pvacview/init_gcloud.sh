LAB="griffith-lab"
REPO="shiny-apps"
IMAGE_NAME="pvacview"
TAG_NAME="latest" # do not change this unless you have looked at how your service fetches the tag
REGION="us-central1"
SERVICE_NAME="pvacview"
SERVICE_TAG=$1 # format like: v1-0-0

# check if Docker is running
if ! docker info >/dev/null 2>&1; then
  echo "Docker does not appear to be running. Please start Docker."
  # Exit with a non-zero status to indicate that Docker isn't running
  exit 1
fi

# build the image 
docker build . --tag $REGION-docker.pkg.dev/$LAB/$REPO/$IMAGE_NAME:$TAG_NAME

# deploy: https://cloud.google.com/sdk/gcloud/reference/run/deploy
gcloud --project $LAB run deploy $SERVICE_NAME \
  --image $REGION-docker.pkg.dev/$LAB/$REPO/$IMAGE_NAME:$TAG_NAME \
  --port=3333 \
  --min-instances=0 \
  --max-instances=100 \
  --memory=512Mi \
  --cpu=1 \
  --timeout=300 \
  --allow-unauthenticated \
  --update-env-vars KEY=VALUE --revision-suffix $SERVICE_TAG
