export default {
    template: `
        <div>
            <nav class="navbar navbar-expand-lg navbar-dark bg-dark">
                <a class="navbar-brand px-2" href="#" @click="changeCurrentDirectory">pVACcompare</a>
                <button class="navbar-toggler" type="button" data-bs-toggle="collapse" data-bs-target="#navbarNav" aria-controls="navbarNav" aria-expanded="false" aria-label="Toggle navigation">
                    <span class="navbar-toggler-icon"></span>
                </button>
                <div class="collapse navbar-collapse" id="navbarNav">
                    <ul class="navbar-nav">
                        <li class="nav-item" v-for="item in comparisonItems" :key="item.id" :class="{ active: item.id === currentPageId }">
                            <a class="nav-link" href="#" @click="changeCurrentPage(item.id)">{{ item.file }}</a>
                        </li>
                    </ul>
                </div>
                <div class="dropdown ms-auto px-2" ref="dropdown">
                    <button class="btn btn-secondary dropdown-toggle" type="button" @click="toggleDropdown">
                        {{ currentClass === 1 ? 'MHC Class I' : 'MHC Class II' }}
                    </button>
                    <ul class="dropdown-menu" :class="{ show: isDropdownOpen }">
                        <li v-if="hasMhcClassI"><a class="dropdown-item" href="#" @click="changeClass(1)">MHC Class I</a></li>
                        <li v-if="hasMhcClassII"><a class="dropdown-item" href="#" @click="changeClass(2)">MHC Class II</a></li>
                    </ul>
                </div>
            </nav>
        </div>
    `,
    props: [
        'currentPageId',
        'currentClass',
        'comparisonItems',
        'hasMhcClassI',
        'hasMhcClassII'
    ],

    data() {
        return {
            isDropdownOpen: false,
        };
    },

    methods: {
        changeCurrentDirectory() {
            this.$emit('change-directory');
        },

        changeCurrentPage(id) {
            this.$emit('change-page', id);
        },

        toggleDropdown() {
            this.isDropdownOpen = !this.isDropdownOpen;
        },

        changeClass(classId) {
            this.isDropdownOpen = false;
            this.$emit('change-class', classId);
        },

        closeDropdownIfOutside(event) {
            const dropdown = this.$refs.dropdown;
            if (dropdown && !dropdown.contains(event.target) && this.isDropdownOpen) {
                this.isDropdownOpen = false;
            }
        }
    },
    mounted() {
        document.addEventListener('click', this.closeDropdownIfOutside);
    },

    beforeUnmount() {
        document.removeEventListener('click', this.closeDropdownIfOutside);
    }
};