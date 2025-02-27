export default {
    template: `
        <div v-if="hasUniqueEntries" class="mx-5 mt-5">
            <h1>Entries</h1>
            <div class="mx-4">
                <div class="mt-3" v-for="(entries, section) in unique_entries" :key="section">
                    <template v-if="entries && Object.keys(entries).length">
                        <h4 v-if="entries"><u>{{ section }}</u></h4>
                        <div class="q-pa-md">
                            <q-virtual-scroll
                                type="table"
                                style="max-height: 70vh"
                                :virtual-scroll-item-size="48"
                                :virtual-scroll-sticky-size-start="48"
                                :items="formattedEntries(entries)"
                            >
                                <template v-slot:before>
                                    <thead class="thead-sticky text-left">
                                        <tr>
                                            <th>
                                                ID
                                                <q-tooltip anchor="top middle" self="top middle">
                                                    {{ idFormat }}
                                                </q-tooltip>
                                            </th>
                                            <th v-if="!Array.isArray(entries)">
                                                Number of Hits    
                                            </th>
                                        </tr>
                                    </thead>
                                </template>
                                <template v-slot="{ item: row }">
                                    <tr class="table-item">
                                        <template v-if="row.value">
                                            <td colspan="2">{{ row.value }}</td>
                                        </template>
                                        <template v-else-if="row.hits">
                                            <td>{{ row.id }}</td>
                                            <td>{{ row.hits }}</td>
                                        </template>
                                        <template v-else>
                                            <td colspan="2">No data available</td>
                                        </template>
                                    </tr>
                                </template>
                            </q-virtual-scroll>
                        </div>
                    </template>
                </div>
            </div>
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
        'referenceMatchesData'
    ],
    data() {
        return {
            unique_entries: {},
            idFormat: ""
        };
    },

    computed: {
        hasUniqueEntries() {
            return Object.values(this.unique_entries).some(
                (entries) => Array.isArray(entries) ? entries.length > 0 : Object.keys(entries).length > 0
            );
        }
    },

    watch: {
        comparisonItems: function(items) {
            if (items.length > 0) {
                this.loadEntries();
            }
        },
        currentPageId(newPageId, oldPageId) {
            if (newPageId !== oldPageId) {
                this.loadEntries();
            }
        },
        hasUniqueEntries: function(value) {
            this.$emit('has-unique-entries-changed', value);
        }
    },

    methods: {
        loadEntries() {
            const key = this.currentComparison.key;
            if (!this[key] || !this[key].entries) {
                return null;
            }

            this.idFormat = this[key].id_format;
            const unique_entries = this[key].entries;
            if (key === 'inputYmlData' || key === 'jsonInputData') {
                this.unique_entries = {};
            }
            this.unique_entries = unique_entries;
        },

        formattedEntries(entries) {
            if (typeof entries === 'object' && !Array.isArray(entries)) {
                return Object.entries(entries).map(([id, hits]) => ({ id, hits }));
            } else if (Array.isArray(entries)) {
                return entries.map((value, index) => ({ id: index + 1, value }));
            }
            return entries;
        }
    }
};