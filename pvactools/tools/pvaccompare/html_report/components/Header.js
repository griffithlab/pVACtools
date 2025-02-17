export default {
    template: `
        <div class="mx-5 mt-4">
            <p><strong>Selected Directory:</strong> {{ currentDirectory }}</p>
            <h1><u>{{ currentComparison.name }}</u></h1>
            <div class="mx-3">
                <ul class="list-unstyled">
                    <li><strong>File 1:</strong> {{ inputFile1 }}</li>
                    <li><strong>File 2:</strong> {{ inputFile2 }}</li>
                </ul>
            </div>
        </div>
    `,
    props: [
        'currentDirectory',
        'currentPageId',
        'currentComparison',
        'comparisonItems',
        'aggregatedData',
        'unaggregatedData',
        'inputYmlData',
        'jsonInputData',
        'referenceMatchesData'
    ],
    computed: {
        inputFile1() {
            const item = this.comparisonItems[this.currentPageId];
            return item && this[item.key] ? this[item.key].input_file1 : "Unavailable";
        },

        inputFile2() {
            const item = this.comparisonItems[this.currentPageId];
            return item && this[item.key] ? this[item.key].input_file2 : "Unavailable";
        }
    }
};