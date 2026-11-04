// Screen recording of one slide stepped through its clicks, from the running dev server (localhost:3030).
// Usage: node tools/record_slide.mjs SLIDE CLICKS SECONDS_PER_STEP OUTDIR
import { chromium } from "playwright-chromium";

const [slide, clicks, dwell, outdir, last] = process.argv.slice(2);
const browser = await chromium.launch();
const ctx = await browser.newContext({ viewport: { width: 1280, height: 720 }, recordVideo: { dir: outdir, size: { width: 1280, height: 720 } } });
const page = await ctx.newPage();
await page.goto(`http://localhost:3030/${slide}?embedded=true`, { waitUntil: "networkidle" });
await page.waitForTimeout(Number(dwell) * 1000);
for (let k = 0; k < Number(clicks); k++) {
  await page.keyboard.press("ArrowRight");
  await page.waitForTimeout(Number(dwell) * 1000);
}
await page.waitForTimeout(Number(last || 0) * 1000);
await ctx.close();
await browser.close();
