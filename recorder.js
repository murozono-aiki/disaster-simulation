const canvas = document.getElementById("canvas");
const canvasStream = canvas.captureStream(30);

const mediaRecorder = new MediaRecorder(canvasStream);

mediaRecorder.start();
setTimeout(() => {mediaRecorder.stop()}, 10000);

mediaRecorder.addEventListener("dataavailable", event => {
    console.log(event);
    const videoBlob = event.data;//new Blob([event.data], { type: event.data.type });
    const dataUrl = window.URL.createObjectURL(videoBlob);
    const anchor = document.getElementById("downloadLink");
    anchor.download = `disaster simulation movie`;
    anchor.href = dataUrl;
});