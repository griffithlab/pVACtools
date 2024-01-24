IMAGE_NAME="pvacview"
TAG_NAME="latest" # do not change this unless you have looked at how your service fetches the tag
REGION="us-central1"
LAB="griffith-lab" # also the project name
REPO="shiny-apps"
SERVICE_NAME="pvacview"
SERVICE_TAG=$2 # format like: v1-0-0

# script
COMMAND=$1 #  init or update

function show_help {
  echo "usage: sh gcloud_service.sh COMMAND <TAG>"
  echo ""
  echo "commands:"
  echo "    init       Creates a new service in Cloud Run."
  echo "    update     updates an existing service in Cloud Run."
  echo ""
  echo "arguments:"
  echo "    <TAG>     Positional argument that tags the revision/release of the service (example: v1-0-0)."
  echo ""
}

function validate {
  # check if Docker is running
  if ! docker info >/dev/null 2>&1; then
    echo "ERROR: Docker does not appear to be running. Please start Docker."
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
}

function build_and_push_image {
  # build the image 
  docker build . --tag $REGION-docker.pkg.dev/$LAB/$REPO/$IMAGE_NAME:$TAG_NAME

  # Check the exit code of the docker build command
  if [ $? -ne 0 ]; then
    echo "Error: Docker build failed"
    exit 1
  fi

  # push the image
  docker push $REGION-docker.pkg.dev/$LAB/$REPO/$IMAGE_NAME:$TAG_NAME

  # Check the exit code of the docker push command
  if [ $? -ne 0 ]; then
    echo "Error: Docker push failed"
    exit 1
  fi
}

function init {
  build_and_push_image

  # deploy: https://cloud.google.com/sdk/gcloud/reference/run/deploy
  gcloud --project $LAB run deploy $SERVICE_NAME \
    --image $REGION-docker.pkg.dev/$LAB/$REPO/$IMAGE_NAME:$TAG_NAME \
    --port=3333 \
    --min-instances=0 \
    --max-instances=100 \
    --memory=2Gi \
    --cpu=4 \
    --timeout=300 \
    --allow-unauthenticated \
    --update-env-vars KEY=VALUE --revision-suffix $SERVICE_TAG
}

function update {
  # check SERVICE_TAG does not already exist
  TAGS=$(gcloud run revisions list --project $LAB --service $SERVICE_NAME --format="value(REVISION)")
  if [[ $TAGS =~ $SERVICE_TAG ]]; then
    echo "ERROR: The Service \"$SERVICE_NAME\" already has revision $SERVICE_TAG. Please use a new tag."
    exit 1
  fi

  # set region; updating service will ask for region if this is not set
  if [[ $(gcloud config get-value run/region) != "$REGION" ]]; then
    echo "<< Upadting property [run/region] to $REGION >>"
    gcloud config set run/region $REGION
  fi

  build_and_push_image

  # update service
  if [[ -z $(gcloud --project $LAB run services list | grep -w $SERVICE_NAME) ]]; then
    echo "<<ERROR: Your service does not exist. You must initialize it first. >>"
  else
    gcloud --project $LAB run services update $SERVICE_NAME --update-env-vars KEY=VALUE --revision-suffix $SERVICE_TAG
  fi
}

if [[ ($COMMAND != "init") && ($COMMAND != "update")]]; then
  show_help
  echo "ERROR: invalid command."
  exit 1
else
  validate
  case $COMMAND in
    "init")
      init
      ;;
    "update")
      update
      ;;
  esac
fi














