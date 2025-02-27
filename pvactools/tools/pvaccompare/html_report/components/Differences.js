export default {
    template: `
        <div v-if="Object.keys(differences).length" class="mx-5 mt-5">
            <h1>Differences</h1>
            <div class="mx-4">
                <div class="mt-3" v-for="(entries, section) in differences" :key="section">
                    <template v-if="entries && Object.keys(entries).length">
                        <h4 v-if="entries"><u>{{ section }}</u></h4>
                        <div v-if="tableNeeded" class="q-pa-md">
                            <q-virtual-scroll
                                type="table"
                                style="max-height: 70vh"
                                :virtual-scroll-item-size="48"
                                :virtual-scroll-sticky-size-start="48"
                                :items="visibleEntries[section] || []"
                                @virtual-scroll="(e) => onVirtualScroll(section, e)"
                            >
                                <template v-slot:before>
                                    <thead class="thead-sticky text-left">
                                        <tr>
                                            <th v-for="field in fields" :key="field">
                                                <q-tooltip v-if="field === 'ID'" anchor="top middle" self="top middle">
                                                    {{ idFormat }}
                                                </q-tooltip>
                                                {{ field }}
                                            </th>
                                        </tr>
                                    </thead>
                                </template>
                                <template v-slot="{ item: row }">
                                    <tr class="table-item" v-for="entry in [row]" :key="entry.ID">
                                        <td v-for="(value, field) in entry" :key="field">
                                            {{ value ? value: "NA" }}
                                        </td>
                                    </tr>
                                </template>
                            </q-virtual-scroll>
                        </div>
                        <div v-else>
                            <ul>
                                <li v-for="(value, field) in entries" :key="field">
                                    <template v-if="Array.isArray(entries)">
                                        {{ value }}
                                    </template>
                                    <template v-else-if="typeof value === 'object' && !Array.isArray(value)">
                                        {{ field }}
                                        <ul>
                                            <li v-for="(subValue, subField) in value" :key="subField">
                                                {{ subField }}: {{ subValue }}
                                            </li>
                                        </ul>
                                    </template>
                                    <template v-else>
                                        {{ field ? field + ': ' + value : value }}
                                    </template>
                                </li>
                            </ul>
                        </div>
                    </template>
                </div>
            </div>
        </div>
        <div v-else-if="!hasUniqueEntries" class="identical-msg mx-5 mt-5">
            <h4>The files are identical</h4>
        </div>
    `,
    props: [
        'currentPageId',
        'currentComparison',
        'comparisonItems',
        'aggregatedData',
        'unaggregatedData',
        'inputYmlData',
        'jsonInputData',
        'referenceMatchesData',
        'hasUniqueEntries'
    ],

    data() {
        return {
            visibleEntries: {},
            loadedSections: {},
            currentSectionIndex: {},
            loadBatchSize: 1000,
            fields: ["Entry #", "ID", "File 1 Value", "File 2 Value", "File 1 Line", "File 2 Line"],
            idFormat: "",
            idNum: 1
        };
    },

    computed: {
        differences() {
            return this.getDifferences();
        },

        tableNeeded() {
            if (this.currentComparison !== "Unavailable") {
                return this.currentComparison.key !== 'inputYmlData' && this.currentComparison.key !== 'jsonInputData';
            }
            return false;
        }
    },

    watch: {
        comparisonItems: function(items) {
            if (items.length > 0) {
                this.initializeData();
            }
        },
        currentPageId(newPageId, oldPageId) {
            if (newPageId !== oldPageId) {
                this.initializeData();
            }
        }
    },

    methods: {
        getDifferences() {
            if (this.currentComparison === "Unavailable") {
                return {};
            }

            const key = this.currentComparison.key;
            if (!this[key] || !this[key].differences) {
                return {};
            }

            const differences = this[key].differences;
            const formattedDifferences = {};

            for (const section in differences) {
                const sectionData = differences[section];

                if (Array.isArray(sectionData)) {
                    formattedDifferences[section] = sectionData;
                } else if (typeof sectionData === 'object') {
                    formattedDifferences[section] = {};

                    for (const nestedKey in sectionData) {
                        const nestedData = sectionData[nestedKey];

                        if (Array.isArray(nestedData)) {
                            formattedDifferences[section][nestedKey] = nestedData;
                        } else if (typeof nestedData === 'object') {
                            formattedDifferences[section][nestedKey] = nestedData;
                        } else {
                            formattedDifferences[section][nestedKey] = nestedData;
                        }
                    }
                }
            }

            return formattedDifferences;
        },

        initializeData() {
            if (this.currentComparison === "Unavailable") {
                return;
            }

            const key = this.currentComparison.key;
            if (this[key] && this[key].differences) {
                this.idFormat = this[key].id_format;
                this.visibleEntries = {};
                this.loadedSections = {};
                this.currentSectionIndex = {};

                for (const section in this[key].differences) {
                    this.visibleEntries[section] = [];
                    this.loadedSections[section] = 0;
                    this.currentSectionIndex[section] = 1;

                    this.loadSection(section);
                }
            }
        },

        loadSection(section) {
            const key = this.currentComparison.key;
            const sectionData = this[key]?.differences?.[section];
            if (!sectionData || this.currentSectionIndex[section] > sectionData.num_sections) {
                return;
            }

            const sectionKey = `section${this.currentSectionIndex[section]}`;
            const newEntries = sectionData[sectionKey] || [];
            this.visibleEntries[section] = [...this.visibleEntries[section], ...newEntries];
            this.currentSectionIndex[section]++;
        },

        onVirtualScroll(section, e) {
            const bottomReached = e.index >= this.visibleEntries[section]?.length - 1;
            if (bottomReached) {
                this.loadSection(section);
            }
        },
    }
};