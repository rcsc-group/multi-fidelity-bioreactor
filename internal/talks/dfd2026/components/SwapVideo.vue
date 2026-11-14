<!-- Loop video that hands over to a release video with no jump: the release file starts with the same frames as the loop,
     so on `on` it starts at the loop's current time and plays on through the release. -->
<script setup>
import { ref, watch, onMounted } from 'vue'
const props = defineProps({ loop: String, rel: String, on: Boolean })
const a = ref(null)
const b = ref(null)
const showB = ref(false)
function swap(v) {
  if (!a.value || !b.value) return
  if (v) {
    b.value.currentTime = a.value.currentTime
    b.value.play()
    showB.value = true
  } else {
    showB.value = false
    b.value.pause()
    a.value.play()
  }
}
watch(() => props.on, swap)
onMounted(() => { if (props.on) { b.value.currentTime = 0; b.value.play(); showB.value = true } })
</script>

<template>
  <div class="swap">
    <video ref="a" :src="loop" autoplay loop muted playsinline :style="{ visibility: showB ? 'hidden' : 'visible' }" />
    <video ref="b" :src="rel" muted playsinline preload="auto" class="top" :style="{ visibility: showB ? 'visible' : 'hidden' }" />
  </div>
</template>

<style scoped>
.swap { position: relative; width: 100%; }
.swap video { width: 100%; display: block; }
.swap .top { position: absolute; inset: 0; }
</style>
