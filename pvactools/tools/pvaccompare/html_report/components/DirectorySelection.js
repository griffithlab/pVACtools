export default {
    template: `
        <div>
            <div class="bg-dark text-white text-left py-4 px-2 shadow-sm">
                <h1 class="display-4 mb-0"><strong>pVACcompare</strong></h1>
            </div>
            <div class="mx-5 mt-5">
                <h2>Select Results Folder</h2>
                <p class="confirm-access-msg"><i>You will be asked to confirm directory access. This is a browser security feature.</i></p>
                <div class="input-group mb-3 mt-4">
                    <input
                        ref="directoryInput"
                        type="file"
                        webkitdirectory
                        directory
                        class="form-control custom-border"
                        @change="onDirectorySelected"
                    />
                </div>

                <div v-if="selectedDirectory" class="mt-3">
                    <q-card flat>
                        <q-card-section>
                            <div class="text-h6">Selected Directory:</div>
                            <div class="text-subtitle2">{{ selectedDirectory }}</div>
                        </q-card-section>
                    </q-card>
                    <q-card flat>
                        <q-card-section v-if="filesInDirectory.length > 0">
                            <div class="text-h6">Files Found in Directory:</div>
                            <ul>
                                <li v-for="file in filesInDirectory" :key="file.id">{{ file.name }}</li>
                            </ul>
                        </q-card-section>
                    </q-card>
                </div>

                <q-btn
                    v-if="selectedDirectory"
                    class="q-mt-md"
                    label="Confirm and Load Files"
                    color="primary"
                    @click="confirmSelection"
                />
            </div>
        </div>
    `,
    data() {
        return {
            filesInDirectory: [],
            validFilesInDirectory: [],
            selectedDirectory: null,
        };
    },
    methods: {
        onDirectorySelected(event) {
            const files = event.target.files;

            if (files.length > 0) {
                const path = files[0].webkitRelativePath.split('/')[0];
                this.selectedDirectory = path;

                const fileKeyMap = {
                    'yml_input_data.json': 'inputYmlData',
                    'unaggregated_data.json': 'unaggregatedData',
                    'aggregated_data.json': 'aggregatedData',
                    'reference_matches_data.json': 'referenceMatchesData',
                    'json_input_data.json': 'jsonInputData',
                };

                Array.from(files).forEach((file, index) => {
                    const key = fileKeyMap[file.name] || null;
                    const reader = new FileReader();
                    reader.onload = (e) => {
                        this.filesInDirectory.push({
                            id: index,
                            name: file.name,
                            content: e.target.result,
                            key,
                        });
                    };
                    reader.readAsText(file);
                });
            } else {
                this.$q.notify({
                    type: 'negative',
                    message: `Error: The selected directory is empty`,
                });
                this.selectedDirectory = null;
                this.filesInDirectory = [];
            }
        },

        confirmSelection() {
            if (!this.selectedDirectory) {
                this.$q.notify({
                    type: 'negative',
                    message: 'Error: Please select a directory',
                });
                return;
            }

            const validFileNames = [
                'aggregated_data.json',
                'json_input_data.json',
                'reference_matches_data.json',
                'unaggregated_data.json',
                'yml_input_data.json'
            ]

            this.validFilesInDirectory = this.filesInDirectory.filter(file => validFileNames.includes(file.name));
            if (this.validFilesInDirectory.length === 0) {
                this.$q.notify({
                    type: 'negative',
                    message: `Error: No valid files were found in the selected directory`,
                });
                return;
            }
            const mhcClassI = [];
            const mhcClassII = [];

            for (const file of this.validFilesInDirectory) {
                if (file.content) {
                    const jsonData = JSON.parse(file.content);
                    if (jsonData.mhc_class === "1") {
                        if (mhcClassI.some(existingFile => existingFile.name === file.name)) {
                            this.$q.notify({
                                type: 'negative',
                                message: `Error: Duplicate file detected in MHC Class I - ${file.name}`,
                            });
                            return;
                        }
                        mhcClassI.push(file);
                    } else if (jsonData.mhc_class === "2") {
                        if (mhcClassII.some(existingFile => existingFile.name === file.name)) {
                            this.$q.notify({
                                type: 'negative',
                                message: `Error: Duplicate file detected in MHC Class II - ${file.name}`,
                            });
                            return;
                        }
                        mhcClassII.push(file);
                    } else {
                        console.warning(`File ${file.name} does not contain a valid mhc_class`);
                    }
                } else {
                    console.warn(`File ${file.name} has no content.`);
                }
            }

            this.$emit('directory-selected', this.selectedDirectory, this.validFilesInDirectory, mhcClassI, mhcClassII);
            this.$q.notify({
                type: 'positive',
                message: 'Files loaded successfully',
            });
        },
    },
};