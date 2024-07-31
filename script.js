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

// カメラを作成
const camera = new THREE.PerspectiveCamera(45, width / height, 1, 10000);
// カメラの初期座標を設定（X座標:0, Y座標:0, Z座標:0）
camera.position.set(0, 0, 1500);

// 箱を作成
const geometry = new THREE.BoxGeometry(500, 500, 500);
const material = new THREE.MeshStandardMaterial({color: 0x0000FF});
const box = new THREE.Mesh(geometry, material);
//scene.add(box);

// 平行光源
const light = new THREE.DirectionalLight(0xFFFFFF);
light.intensity = 2; // 光の強さを倍に
light.position.set(1, 1, 1); // ライトの方向
// シーンに追加
scene.add(light);


let rot = 0; // 角度
let mouseX = 0; // マウス座標

// マウス座標はマウスが動いた時のみ取得できる
document.addEventListener("mousemove", (event) => {
  mouseX = event.pageX;
});


// 初回実行
tick();

function tick() {

    // 箱を回転させる
    box.rotation.x += 0.01;
    box.rotation.y += 0.01;

    // レンダリング
    renderer.render(scene, camera);


    // マウスの位置に応じて角度を設定
    // マウスのX座標がステージの幅の何%の位置にあるか調べてそれを360度で乗算する
    const targetRot = (mouseX / window.innerWidth) * 360;
    // イージングの公式を用いて滑らかにする
    // 値 += (目標値 - 現在の値) * 減速値
    rot += (targetRot - rot) * 0.02;

    // ラジアンに変換する
    const radian = rot * Math.PI / 180;
    // 角度に応じてカメラの位置を設定
    camera.position.x = 1000 * Math.sin(radian);
    camera.position.z = 1000 * Math.cos(radian);
    // 原点方向を見つめる
    camera.lookAt(new THREE.Vector3(0, 0, 0));

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


function finishRecord() {
    const event = new Event("finishRender");
    document.getElementById("canvas").dispatchEvent(event);
}