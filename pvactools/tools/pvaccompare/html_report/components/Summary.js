export default {
    template: `
        <div class="mx-5 mt-4">
            <div v-if="displaySummary" class="d-flex align-items-center">
                <h3>Summary</h3>
                <button class="btn btn-primary btn-sm m-2" @click.prevent="toggleHideSummary">
                    {{ hide_summary ? 'Show' : 'Hide' }}
                </button>
            </div>
            <div v-if="!hide_summary" class="mx-4">
                <div v-for="(entries, section) in summary" :key="section">
                    <template v-if="Object.keys(entries).length > 0">
                        <div class="summary-notes" v-if="section === 'Notes'">
                            <ul>
                                <li v-for="(value, index) in entries" :key="index">{{ value }}</li>
                            </ul>
                        </div>
                        <div v-else>
                            <h4><u>{{ section }}</u></h4>
                            <ul>
                                <li v-for="(value, field) in entries" :key="field">{{ field }}: {{ value }}</li>
                            </ul>
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
            hide_summary: false,
            summary: {}
        };
    },

    computed: {
        displaySummary() {
            return this.currentComparison.key !== 'inputYmlData' && this.currentComparison.key !== 'jsonInputData';
        }
    },

    watch: {
        comparisonItems: function(items) {
            if (items.length > 0) {
                this.loadSummary();
            }
        },
        currentPageId(newPageId, oldPageId) {
            if (newPageId !== oldPageId) {
                this.loadSummary();
            }
        }
    },

    methods: {
        toggleHideSummary() {
            this.hide_summary = !this.hide_summary;
        },

        loadSummary() {
            const key = this.currentComparison.key;
            if (!this[key] || !this[key].summary) {
                return null;
            }

            const summary = this[key].summary;
            if (key === 'inputYmlData') {
                this.summary = {};
            } else if (key == 'jsonInputData') {
                this.summary = {};
            }
            this.summary = summary;
        }
    }
};