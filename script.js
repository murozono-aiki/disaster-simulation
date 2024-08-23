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

// マウス座標はマウスが動いた時のみ取得できる
document.addEventListener("mousemove", onMousemove_camera);

// 非同期処理で待機するのでasync function宣言とする
async function init() {
    // GLTF形式のモデルデータを読み込む
    const loader = new GLTFLoader();
    // GLTFファイルのパスを指定
    const gltf = await loader.loadAsync('data/Juso-data.glb');
    // 読み込み後に3D空間に追加
    const model = gltf.scene;
    model.position.set(0, -10, 0);
    scene.add(model);
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


/* シミュレーション部分 */

const h = 0.012; //影響半径
const particleMass = 0.0002; //粒子の質量

const densityCoef = particleMass * 315 / (64 * Math.PI * Math.pow(h,9)); //密度計算で使うヤツ
const pressureCoef = particleMass * 45 / (Math.PI * Math.pow(h,6)); //圧力計算で使うヤツ

const pressureStiffness = 200; //圧力係数
const restDensity = 1000; //静止密度
 
/**
 * 粒子の密度計算
 * @param {{position:THREE.Vector3, velocity:THREE.Vector3, force:THREE.Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function CalcDensity(particles) {
    const h2 = h*h; //事前にhの二乗を計算しておく
    for (let i = 0; i < particles.length; i++) { //一つづつ粒子の密度を計算
        let nowParticle = particles[i]; //今回計算する粒子
        let sum = 0; //足し合わせる変数
        for (let j = 0; j < particles.length; j++) { //他の粒子全てについて
            if(i == j){continue;} //自分自身だったらスキップ
            let nearParticle = particles[j];
            
            let diff = nearParticle.position.clone().sub(nowParticle.position.clone()); //粒子距離
            const r2 = diff.clone().dot(diff.clone()); //粒子距離の２乗
 
            //粒子距離がhより小さい場合だけ計算する
            if ( r2 < h2 ) {
                const c = h2 - r2;
                sum += Math.pow(c,3); //(h2-r2)の３乗
            }
        }
 
        nowParticle.density = sum * densityCoef; //密度が求まった
    }
}

/**
 * 粒子の圧力計算
 * @param {{position:THREE.Vector3, velocity:THREE.Vector3, force:THREE.Vector3, density:number, pressure:number}[]} particles - 粒子のリスト
 */
function CalcPressure(particles) {
    for(let i = 0; i < particles.length; i++) { //一つづつ粒子の圧力を計算
        particles[i].pressure = pressureStiffness * ( particles[i].density - restDensity );
    }
}


function tick() {

    

    // レンダリング
    renderer.render(scene, camera);

    requestAnimationFrame(tick);
}

// 初回実行
tick();