LAB="griffith-lab"
REPO="shiny-apps"
IMAGE_NAME="pvacview"
TAG_NAME="latest"
SERVICE_NAME="pvacview-latest"
REGION="us-central1"

# check if Docker is running
if ! docker info >/dev/null 2>&1; then
  echo "Docker does not appear to be running. Please start Docker."
  # Exit with a non-zero status to indicate that Docker isn't running
  exit 1
fi

# check if you are logged in
if [[ -z $(gcloud auth list --filter=status:ACTIVE --format="value(account)") ]]; then
  echo "<< Please log in.>>"
  gcloud auth login
else
  echo "<< [Alread logged in]: Welcome, $(gcloud auth list --filter=status:ACTIVE --format="value(account)") >>"
fi

# check if user is allowed to push to Artifact Registry
if [[ -z $(cat ~/.docker/config.json | grep "\"$REGION-docker.pkg.dev\": \"gcloud\"") ]]; then
  echo "<< Attempting to Authenticate before pushing to Artifact Registry. >>"
  gcloud auth configure-docker $REGION-docker.pkg.dev
else
  echo "<< Allowed to push to Artifact Registry  >>"
fi

# build the image 
docker build . --tag $REGION-docker.pkg.dev/$LAB/$REPO/$IMAGE_NAME:$TAG_NAME

# push the image
docker push $REGION-docker.pkg.dev/$LAB/$REPO/$IMAGE_NAME:$TAG_NAME

# set region; updating service will ask for region if this is not set
if [[ $(gcloud config get-value run/region) != "$REGION" ]]; then
  echo "<< Upadting property [run/region] to $REGION >>"
  gcloud config set run/region $REGION
fi

# update service
if [[ -z $(gcloud --project $LAB run services list | grep -w $SERVICE_NAME) ]]; then
  echo "<< Your service does not exist. You must initialize it first. >>"
else
  gcloud --project $LAB run services update $SERVICE_NAME --update-env-vars KEY=VALUE
fi