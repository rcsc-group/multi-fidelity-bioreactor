<template>
  <div v-if="$nav.currentPage > 1" class="slnum-global">
    {{ $nav.currentPage }}
  </div>
</template>

<style scoped>
.slnum-global {
  position: absolute;
  right: 18px;
  bottom: 10px;
  font-family: "Helvetica Neue", "Inter", sans-serif;
  font-weight: 200;
  font-size: 20px;
  color: var(--ink, #1f2933);
  z-index: 50;
}
</style>
