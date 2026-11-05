<template>
  <div v-if="$nav.currentPage > 1" class="slnum-global">
    {{ $nav.currentPage }}
  </div>
</template>

<style scoped>
.slnum-global {
  position: absolute;
  right: 14px;
  bottom: 8px;
  font-size: 12px;
  color: var(--mute, #7a8797);
  z-index: 50;
}
</style>
