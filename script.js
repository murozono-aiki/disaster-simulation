import * as THREE from "three";
import { GLTFLoader } from "three/addons/loaders/GLTFLoader.js"

const width = 1280;
const height = 720;

// レンダラーを作成
const renderer = new THREE.WebGLRenderer({
    canvas: document.getElementById("canvas")
});
renderer.setSize(width, height);
renderer.setPixelRatio(window.devicePixelRatio);

// シーンを作成
const scene = new THREE.Scene();
scene.background = new THREE.Color(0xFFFFFF);

// カメラを作成
const camera = new THREE.PerspectiveCamera(45, width / height, 1, 10000);
camera.position.set(300, 100, 1000);

// 球を作成
const geometry = new THREE.SphereGeometry(1);
const material = new THREE.MeshStandardMaterial({color: 0x0000FF});
const sphere = new THREE.Mesh(geometry, material);
sphere.position.y = 100;
scene.add(sphere);

// 平行光源
const light = new THREE.DirectionalLight(0xFFFFFF);
light.intensity = 2; // 光の強さを倍に
light.position.set(1, 1, 1); // ライトの方向
// 環境光
const ambientLight = new THREE.AmbientLight(0xbbbbbb);
scene.add(ambientLight);
// シーンに追加
scene.add(light);


let rot = 0; // 角度
let mouseX = 0; // マウス座標
let mouseY = 0; // マウス座標

function onMousemove_camera(event) {
    const currentMouseX = event.pageX;
    const currentMouseY = event.pageY;

    if (event.buttons % 2 == 1) {  // 左クリック時
        const cameraWorldDirection = new THREE.Vector3();
        camera.getWorldDirection(cameraWorldDirection);
        const cameraDirection = cameraWorldDirection.normalize();
        const cameraPositionVectorY = new THREE.Vector3(0, 1, 0);
        const cameraPositionVectorX = cameraDirection.clone().cross(cameraPositionVectorY.clone()).normalize();
        if (!event.ctrlKey) {
            const cameraMoveX = cameraPositionVectorX.clone().multiplyScalar(currentMouseX - mouseX);
            camera.position.x -= cameraMoveX.x;
            camera.position.y -= cameraMoveX.y;
            camera.position.z -= cameraMoveX.z;
            if (!event.shiftKey) {
                const cameraMoveY = cameraPositionVectorY.clone().multiplyScalar(currentMouseY - mouseY);
                camera.position.x += cameraMoveY.x;
                camera.position.y += cameraMoveY.y;
                camera.position.z += cameraMoveY.z;
            } else {
                const cameraMoveZ = cameraDirection.clone().multiplyScalar(currentMouseY - mouseY);
                camera.position.x += cameraMoveZ.x;
                camera.position.y += cameraMoveZ.y;
                camera.position.z += cameraMoveZ.z;
            }
        } else {
            if (Math.abs(currentMouseY - mouseY) > Math.abs(currentMouseX - mouseX)) {
                camera.rotateOnAxis(new THREE.Vector3(1, 0, 0), (currentMouseY - mouseY) / 80);
            } else {
                camera.rotateOnWorldAxis (new THREE.Vector3(0, 1, 0), (currentMouseX - mouseX) / 80);
            }
        }
    }

    mouseX = currentMouseX;
    mouseY = currentMouseY;
}

function onMousedown_camera(event) {}

// マウス座標はマウスが動いた時のみ取得できる
document.addEventListener("mousemove", onMousemove_camera);


// 初回実行
tick();

function tick() {

    // 箱を回転させる
    //box.rotation.x += 0.01;
    //box.rotation.y += 0.01;

    // マウスの位置に応じて角度を設定
    // マウスのX座標がステージの幅の何%の位置にあるか調べてそれを360度で乗算する
    //const targetRot = (mouseX / window.innerWidth) * 360;
    // イージングの公式を用いて滑らかにする
    // 値 += (目標値 - 現在の値) * 減速値
    //rot += (targetRot - rot) * 0.02;

    // ラジアンに変換する
    //const radian = rot * Math.PI / 180;
    // 角度に応じてカメラの位置を設定
    //camera.position.x = 1000 * Math.sin(radian);
    //camera.position.z = 1000 * Math.cos(radian);
    // 原点方向を見つめる
    //camera.lookAt(new THREE.Vector3(0, 0, 0));

    // レンダリング
    renderer.render(scene, camera);

    requestAnimationFrame(tick);
}

// 非同期処理で待機するのでasync function宣言とする
async function init() {
    // ･･･省略
  
    // GLTF形式のモデルデータを読み込む
    const loader = new GLTFLoader();
    // GLTFファイルのパスを指定
    const gltf = await loader.loadAsync('data/Juso-data.glb');
    // 読み込み後に3D空間に追加
    const model = gltf.scene;
    model.position.set(0, -10, 0);
    scene.add(model);
    
    // ･･･省略
}
init();


function startRecord() {
    const event = new Event("startRender");
    document.getElementById("canvas").dispatchEvent(event);
}
function finishRecord() {
    const event = new Event("finishRender");
    document.getElementById("canvas").dispatchEvent(event);
}